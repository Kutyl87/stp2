%% Parametry regulatora DM
N = 22;
Nu = 2;
D = 90; 
%% Parametry równania transmitancji
T0 = 5;
T1= 1.78;
T2 = 5.13;
Tp = 0.5;
T_0base = 5;
K_base = 4.7;
kk=1000;
Gs = tf(K_base,[T1*T2, T1+T2, 1],'IODelay',T0);
Gz = c2d(Gs,Tp,"zoh");
s = step(Gz,0:Tp:100);
%% Inicjalizacja wektorów
u(1:kk)=0; y(1:kk)=0;
yzad(1:kk)=0; yzad(abs(-2- T0 *(1/Tp))+1:kk)=1;
e(1:kk)=0;
T0_prop = [1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2];
K_prop_dmc = [];
lambda = 1;
Gamma = eye(N,N);
k_start = 1000;
Alpha = eye(Nu,Nu) * lambda;
%% Główna pętla do poszukiwania granicy stabilności
for x= 1 : size(T0_prop,2)
    T0 = T_0base * T0_prop(x);
    z = 0.01;
    left_max = 0;
    right_max = 0;
    left_min = 0;
    right_min = 0;
    while ((left_max*1.001>= right_max & left_min<= right_min*1.001)  | k_start >800)
        k_start = 900;
        K_t = K_base*z;
        Gs = tf(K_t,[T1*T2, T1+T2, 1],'IODelay',T0);
        Gz = c2d(Gs,Tp,"zoh");
        a1 = Gz.Denominator{1}(2);
        a0 = Gz.Denominator{1}(3);
        b1 = Gz.Numerator{1}(2);
        b0 = Gz.Numerator{1}(3);
        M = zeros(N,Nu);
        Mp = zeros(N,D-1);
        kk=1000; 
        u(1:kk)=0; y(1:kk)=0;
        yzad(1:kk)=0; yzad(abs(-2- T0 *(1/Tp))+1:kk)=1;
        e(1:kk)=0;
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
     if (y(k) >= yzad(k)) && k_start == 900
        k_start = k;
    end
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
    z = z+0.001;
    left_max = max(y(k_start+1:min(k_start+100,1000)));
    left_min = min(y(k_start+1:min(k_start+100,1000)));
    right_max = max(y(800:1000));
    right_min = min(y(800:1000));
    end
    K_prop_dmc = [K_prop_dmc, z];
end