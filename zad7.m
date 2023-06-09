K =4.7;
T0 = 5;
T1= 1.78;
T2 = 5.13;
Tp = 0.5;
Gs = tf(K,[T1*T2, T1+T2, 1],'IODelay',T0);
Gz = c2d(Gs,Tp,"zoh");
N = 90;
Nu = 90;
D = 90;
M = zeros(N,Nu);
Mp = zeros(N,D-1);
s = step(Gz,0:Tp:100);
kk=1000;
u(1:D-1)=0; y(1:D-1)=0;
yzad(1:D-1)=0; yzad(D:kk)=1;
e(1:D-1)=0;
for i= 1:N
    for j = 1:Nu
        if (i-j+1)>0
            M(i,j) = s(i-j+1);
        end
    end
end

% for i = 1:N
%     for j = 1:D-1
%         Mp(i,j) = s(i+j) - s(j);
%     end
% end
lambda = 1;
Gamma = eye(N,N);
Alpha = eye(Nu,Nu) * lambda;

K = inv(M' * Gamma * M + Alpha) * M' * Gamma;



du = zeros(0:12);



a1 = Gz.Denominator{1}(2);
a2 = Gz.Denominator{1}(3);
b11 = Gz.Numerator{1}(2);
b12 = Gz.Numerator{1}(3);
coeffs_b = [zeros(1,10),b11,b12];
coeffs_a = [a1,a2];
counter = 2;
predicted_y = 0;
d= zeros(N);
for k=D:kk
 dUp = [];
 y(k)=b11*u(k-1- T0 *(1/Tp))+b12*u(k-2- T0 *(1/Tp))-a1*y(k-1)-a2*y(k-2);
 e(k)=yzad(k)-y(k);
 Yzadk = yzad(k) *ones(N,1);
 Y0 = [];
 d(k) = y(k) - predicted_y;
 for p=1:N
     N_un = min(p,abs(-2- T0 *(1/Tp)));
     N_y = min(p-1,2);
     y0_pr = 0;
     for j=1:N_un
         y0_pr = y0_pr + coeffs_b(j)*u(k-1);
     end
     for j=N_un+1:abs(-2- T0 *(1/Tp))
         y0_pr = y0_pr + coeffs_b(j)*u(k-1);
     end
     for j=1:N_y
         y0_pr = y0_pr -coeffs_a(j)*Y0(end-j+1);
     end
     for j=N_y+1:2
         y0_pr = y0_pr - coeffs_a(j)*y(k-j+p);
     end
     y0_pr = y0_pr + d(k);
     Y0 = [Y0,y0_pr];
 end
 predicted_y = Y0(1);

 % for i=1:D-1
 %     if (k-i-1) > 0
 %        dUp = [dUp;u(k-i) - u(k-i-1)];
 %     else
 %        dUp = [dUp;u(k-i)]
 %     end
 % end
 dU = K*(Yzadk - Y0);
 u(k)= dU(1) + u(k-1);
 % if mod(k,2*D-1) == 0
 %    yzad(k+1:kk)= (1 - min(yzad(k),1))*counter;
 %    counter=counter+1;
 % end 

end
t = linspace(1,kk,kk);
figure; 
stairs(t,u,'LineWidth',1.5, Color='r');
title('u - sterowanie'); 
xlabel('k - number próbki');
ylabel("Wartość sterowania")
figure; 
stairs(t,y,'LineWidth',1.5); 
hold on;
stairs(t,yzad,'LineWidth',1, 'LineStyle','--');
title('Charakterystyki y,y_{zad}'); 
xlabel('k - number próbki');
ylabel('Wartość')
legend("Wartość na wyjściu y", "Wartość zadana y_{zad}",Location="southeast")