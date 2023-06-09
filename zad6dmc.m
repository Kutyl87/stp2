N = 22
Nu = 2
D = 90
coeff = 0;
kk=1000; 
T0 = 5;
T1= 1.78;
T2 = 5.13;
Tp = 0.5;
T_0base = 5;
K_base = 4.7
kk=1000; 
u(1:kk)=0; y(1:kk)=0;
yzad(1:kk)=0; yzad(abs(-2- T0 *(1/Tp))+1:kk)=1;
e(1:kk)=0;
T0_prop = [1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2];
K_prop = [];
lambda = 1;
Gamma = eye(N,N);
Alpha = eye(Nu,Nu) * lambda;
for x= 1 : size(T0_prop,2)
    T0 = T_0base * T0_prop(x);
    z = 1;
    left_max = 0;
    right_max = 0;
    while left_max>= right_max
        M = zeros(N,Nu);
        Mp = zeros(N,D-1);
        kk=1000; 
        u(1:kk)=0; y(1:kk)=0;
        yzad(1:kk)=0; yzad(abs(-2- T0 *(1/Tp))+1:kk)=1;
        e(1:kk)=0;
        K_t = K_base*z;
        Gs = tf(K_t,[T1*T2, T1+T2, 1],'IODelay',T0);
        Gz = c2d(Gs,Tp,"zoh");
        a1 = Gz.Denominator{1}(2);
        a0 = Gz.Denominator{1}(3);
        b1 = Gz.Numerator{1}(2);
        b0 = Gz.Numerator{1}(3);
        s = step(Gz,0:Tp:100);
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
        K = inv(M' * Gamma * M + Alpha) * M' * Gamma;
    for k=D:kk
     dUp = [];
     y(k)=b1*u(k-1- T0 *(1/Tp))+b0*u(k-2- T0 *(1/Tp))-a1*y(k-1)-a0*y(k-2);
     e(k)=yzad(k)-y(k);
     Yzadk = yzad(k) *ones(N,1);
     Yk = y(k) *ones(N,1);
     for i=1:D-1
         if (k-i-1) > 0
            dUp = [dUp;u(k-i) - u(k-i-1)];
         else
            dUp = [dUp;u(k-i)];
         end
     end
    dU = K*(Yzadk - Yk - Mp * dUp);
    u(k)= dU(1) + u(k-1);
    end
    z = z-0.001;
    left_max = max(u(200:400));
    right_max = max(u(800:1000));
    end
    K_prop = [K_prop, z];
end
t = linspace(1,kk,kk)
stairs(t,u,'LineWidth',1.5, Color='r');