%%Parametry regulatora PID
Kp = 0.18;
Ti = 7.95;
Td = 1.94;

%%Wyznaczone wartości r1,r2,r0
r1 = Kp*((Tp/(2*Ti)) -2 *(Td/Tp) -1);
r2 = Kp*Td/Tp;
r0 = Kp*(1+(Tp/(2*Ti)) + (Td/Tp));

%%Parametry równania różnicowego
a1 = Gz.Denominator{1}(2);
a0 = Gz.Denominator{1}(3);
b1 = Gz.Numerator{1}(2);
b0 = Gz.Numerator{1}(3);

%% Inicjalizacja wektorów
kk=1000; 
u(1:kk)=0; y(1:kk)=0;
yzad(1:kk)=0; yzad(abs(-2- T0 *(1/Tp))+1:kk)=1;
e(1:kk)=0;

%% Pętla regulatora
for k=abs(-2- T0 *(1/Tp))+1:kk
 y(k)=b1*u(k-1- T0 *(1/Tp))+b0*u(k-2- T0 *(1/Tp))-a1*y(k-1)-a0*y(k-2);
 e(k)=yzad(k)-y(k);
 u(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
end
%% Przygotowanie wykresów i wizualizacja 
t = linspace(1,kk,kk);
figure
stairs(t,u,'LineWidth',1.5, Color='r');
title('u - sterowanie'); 
xlabel('k - number próbki');
ylabel("Wartość sterowania")
print('zad3poprawu.png','-dpng','-r400')
figure
stairs(t,y,'LineWidth',1.5); 
hold on;
stairs(t,yzad,'LineWidth',1, 'LineStyle','--');
title('Charakterystyki y,y_{zad}'); 
xlabel('k - number próbki');
ylabel('Wartość')
legend("Wartość na wyjściu y", "Wartość zadana y_{zad}",Location="southeast")
print('zad3poprawy.png','-dpng','-r400')