N = 30
Nu = 8
D = 80
M = zeros(N,Nu)
Mp = zeros(N,D-1)
s = step(Gz,0:Tp:100)
kk=300
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

for i = 1:N
    for j = 1:D-1
        Mp(i,j) = s(i+j) - s(j);
    end
end
lambda = 3;
Gamma = eye(N,N)
Alpha = eye(Nu,Nu) * lambda;

K = inv(M' * Gamma * M + Alpha) * M' * Gamma;



du = zeros(0:12)



a1 = Gz.Denominator{1}(2)
a0 = Gz.Denominator{1}(3)
b1 = Gz.Numerator{1}(2)
b0 = Gz.Numerator{1}(3)
for k=D:kk
 dUp = []
 y(k)=b1*u(k-11)+b0*u(k-12)-a1*y(k-1)-a0*y(k-2);
 e(k)=yzad(k)-y(k);
 Yzadk = yzad(k) *ones(N,1);
 Yk = y(k) *ones(N,1);
 for i=1:D-1
     if (k-i-1) > 0
        dUp = [dUp;u(k-i) - u(k-i-1)];
     else
        dUp = [dUp;0]
     end
 end
 dU = K*(Yzadk - Yk - Mp * dUp)
 u(k)= dU(1) + u(k-1);
end;
t = linspace(1,kk,kk)
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