%% Parametry regulatora GPC
N = 22;
Nu = 2;
D = 90;
coeff = 0;

%% Parametry transmitancji
T0 = 5;
T1= 1.78;
T2 = 5.13;
Tp = 0.5;
T_0base = 5;
K_base = 4.7;
kk=1000;
Gs = tf(K_base,[T1*T2, T1+T2, 1],'IODelay',T0);
Gz = c2d(Gs,Tp,"zoh");

%% Wyznaczenie odpowiedzi skokowej
s = step(Gz,0:Tp:100);

%% Inicjalizacja wektorów
u(1:kk)=0; y(1:kk)=0;
yzad(1:kk)=0; yzad(abs(-2- T0 *(1/Tp))+1:kk)=1;
e(1:kk)=0;
T0_prop = [1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2];
K_prop_gpc = [];
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
        a2 = Gz.Denominator{1}(3);
        b11 = Gz.Numerator{1}(2);
        b12 = Gz.Numerator{1}(3);
        coeffs_b = [zeros(1,abs(-2- T0 *(1/Tp))-2),b11,b12];
        coeffs_a = [a1,a2];
        predicted_y = 0;
        d = 0;
        % d= zeros(N,1);
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
        K = inv(M' * Gamma * M + Alpha) * M' * Gamma;
    for k=D:kk
        d = 0;
        dUp = [];
        y(k)=b11*u(k-1- T0 *(1/Tp))+b12*u(k-2- T0 *(1/Tp))-a1*y(k-1)-a2*y(k-2);
        e(k)=yzad(k)-y(k);
        if (y(k) >= yzad(k)) && k_start == 900
            k_start = k;
    end
        Yzadk = yzad(k) *ones(N,1);
        Y0 = [];
        d = y(k) - predicted_y;
        for p=1:N
             N_un = min(p,abs(-2- T0 *(1/Tp)));
             N_y = min(p-1,2);
             y0_pr = 0;
             for j=1:N_un
                 y0_pr = y0_pr + coeffs_b(j)*u(k-1);
             end
             for j=N_un+1:abs(-2- T0 *(1/Tp))
                 y0_pr = y0_pr + coeffs_b(j)*u(k-j+p);
             end
             for j=1:N_y
                 y0_pr = y0_pr -coeffs_a(j)*Y0(end-j+1);
             end
             for j=N_y+1:2
                 y0_pr = y0_pr - coeffs_a(j)*y(k-j+p);
             end
             y0_pr = y0_pr + d;
             Y0 = [Y0,y0_pr];
        end
 predicted_y = Y0(1);
 dU = K*(Yzadk - Y0');
 u(k)= dU(1) + u(k-1);
    end
z = z+0.001;
left_max = max(y(k_start+1:min(k_start+100,1000)));
left_min = min(y(k_start+1:min(k_start+100,1000)));
right_max = max(y(800:1000));
right_min = min(y(800:1000));
    end
    K_prop_gpc = [K_prop_gpc, z];
end