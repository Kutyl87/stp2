K =4.7;
T0 = 5;
T1= 1.78;
T2 = 5.13;
Tp = 0.5;
Gs = tf(K,[T1*T2, T1+T2, 1],'IODelay',T0);
Gz = c2d(Gs,Tp,"zoh");
Kp = 0.18;
Ti = 7.95;
Td = 1.94;
r1 = Kp*((Tp/(2*Ti)) -2 *(Td/Tp) -1);
r2 = Kp*Td/Tp;
r0 = Kp*(1+(Tp/(2*Ti)) + (Td/Tp));
a1 = Gz.Denominator{1}(2);
a0 = Gz.Denominator{1}(3);
b1 = Gz.Numerator{1}(2);
b0 = Gz.Numerator{1}(3);
N = 22;
Nu = 2;
D = 90;
kk=1000; 
u_pid(1:D-1)=0; y_pid(1:D-1)=0;
yzadpid(1:D-1)=0; yzadpid(D:kk)=1;
epid(1:D-1)=0;
counter = 2;
for k=D:kk
    y_pid(k)=b1*u_pid(k-1- T0 *(1/Tp))+b0*u_pid(k-2- T0 *(1/Tp))-a1*y_pid(k-1)-a0*y_pid(k-2);
    epid(k)=yzadpid(k)-y_pid(k);
    u_pid(k)=r2*epid(k-2)+r1*epid(k-1)+r0*epid(k)+u_pid(k-1);
    if mod(k,2*D-1) == 0
        yzadpid(k+1:kk)= (1 - min(yzadpid(k),1))*counter
            counter=counter+1
    end 
end;
t = linspace(1,kk,kk);
M = zeros(N,Nu);
Mp = zeros(N,D-1);
s = step(Gz,0:Tp:100);
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
lambda = 969;
Gamma = eye(N,N)
Alpha = eye(Nu,Nu) * lambda;

K = inv(M' * Gamma * M + Alpha) * M' * Gamma;



du = zeros(0:12)
counter = 2
for k=D:kk
 dUp = []
 y(k)=b1*u(k-1- T0 *(1/Tp))+b0*u(k-2- T0 *(1/Tp))-a1*y(k-1)-a0*y(k-2);
 e(k)=yzad(k)-y(k);
 Yzadk = yzad(k) *ones(N,1);
 Yk = y(k) *ones(N,1);
 for i=1:D-1
     if (k-i-1) > 0
        dUp = [dUp;u(k-i) - u(k-i-1)];
     else
        dUp = [dUp;u(k-i)]
     end
 end
 dU = K*(Yzadk - Yk - Mp * dUp)
 u(k)= dU(1) + u(k-1);
 if mod(k,2*D-1) == 0
    yzad(k+1:kk)= (1 - min(yzad(k),1))*counter
    counter=counter+1
 end 

end;
figure; 
stairs(t,u,'LineWidth',1.5, Color='r');
hold on
stairs(t,u_pid,'LineWidth',1.5, Color='cyan',LineStyle="-.");
title('u - sterowanie'); 
xlabel('k - number próbki');
ylabel("Wartość sterowania")
legend("Regulator DMC", "Regulator PID",Location="southeast")
figure; 
stairs(t,y,'LineWidth',1.5); 
hold on;
stairs(t,y_pid,'LineWidth',1.5); 
stairs(t,yzad,'LineWidth',1, 'LineStyle','--');
title('Charakterystyki y,y_{zad}'); 
xlabel('k - number próbki');
ylabel('Wartość')
legend("y regulatora DMC","y regulatora PID",  "Wartość zadana y_{zad}",Location="southeast")
