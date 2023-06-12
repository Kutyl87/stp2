%% Parametry równania transmitancji
K =4.7;
T0 = 5;
T1= 1.78;
T2 = 5.13;
Tp = 0.5;
Gs = tf(K,[T1*T2, T1+T2, 1],'IODelay',T0);
Gz = c2d(Gs,Tp,"zoh");

%% Parametry regulatora DMC
N = 22;
Nu = 2;
D = 90;

%% Inicjacja macierzy M i Mp
M = zeros(N,Nu);
Mp = zeros(N,D-1);

%% Wyznaczenie współczynników odpowiedzi skokowej
s = step(Gz,0:Tp:100);
kk=1000;
%% Inicjalizacja wektorów
u(1:D-1)=0; y(1:D-1)=0;
yzad(1:D-1)=0; yzad(D:kk)=1;
e(1:D-1)=0;
%% Macierz M
for i= 1:N
    for j = 1:Nu
        if (i-j+1)>0
            M(i,j) = s(i-j+1);
        end
    end
end
%% Macierz Mp
for i = 1:N
    for j = 1:D-1
        Mp(i,j) = s(i+j) - s(j);
    end
end
%% Wyznaczenie K i dobranie parametrów kary
lambda = 969;
Gamma = eye(N,N);
Alpha = eye(Nu,Nu) * lambda;
K = inv(M' * Gamma * M + Alpha) * M' * Gamma;

%% Inicjalizacja wektora du i współczynniki równania różnciowego
du = zeros(0:12);
a1 = Gz.Denominator{1}(2);
a0 = Gz.Denominator{1}(3);
b1 = Gz.Numerator{1}(2);
b0 = Gz.Numerator{1}(3);
counter = 2;
%% Główna pętla regulatora
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
 dU = K*(Yzadk - Yk - Mp * dUp)
 u(k)= dU(1) + u(k-1);
 %% Wprowadzenie zmian wartości zadanej
 if mod(k,2*D-1) == 0
    yzad(k+1:kk)= (1 - min(yzad(k),1))*counter;
    counter=counter+1;
 end 

end

%% Wizualizacja
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