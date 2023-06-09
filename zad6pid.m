coeff = 0;
kk=1000; 
u(1:kk)=0; y(1:kk)=0;
yzad(1:kk)=0; yzad(13:kk)=1;
e(1:kk)=0;
T0 = 5;
T1= 1.78;
T2 = 5.13;
Tp = 0.5;
Kp = 0.18;
Ti = 7.95;
Td = 1.94;
T_0base = 5;
K_base = 4.7
T0_prop = [1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2];
K_prop = [];
for j= 1 : size(T0_prop,2)
    T0 = T_0base * T0_prop(j);
    r1 = Kp*((Tp/(2*Ti)) -2 *(Td/Tp) -1);
    r2 = Kp*Td/Tp;
    r0 = Kp*(1+(Tp/(2*Ti)) + (Td/Tp));
    i = 1;
    left_max = 0;
    right_max = 0;
    while left_max>= right_max
        u(1:kk)=0; y(1:kk)=0;
        yzad(1:kk)=0; yzad(abs(-2- T0 *(1/Tp))+1:kk)=1;
        e(1:kk)=0;
        K = K_base*i;
        Gs = tf(K,[T1*T2, T1+T2, 1],'IODelay',T0);
        Gz = c2d(Gs,Tp,"zoh");
        a1 = Gz.Denominator{1}(2);
        a0 = Gz.Denominator{1}(3);
        b1 = Gz.Numerator{1}(2);
        b0 = Gz.Numerator{1}(3);
        u(1:kk)=0; y(1:kk)=0;
        yzad(1:kk)=0; yzad(abs(-2- T0 *(1/Tp))+1:kk)=1;
        e(1:kk)=0;
        for k=abs(-2- T0 *(1/Tp))+1:kk
         y(k)=b1*u(k-1- T0 *(1/Tp))+b0*u(k-2- T0 *(1/Tp))-a1*y(k-1)-a0*y(k-2);
         e(k)=yzad(k)-y(k);
         u(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
        end;
        i = i+0.001;
        left_max = max(u(200:400));
        right_max = max(u(800:1000));
    end
    K_prop = [K_prop, i];
end