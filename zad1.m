K =20
T0 = 5
T1= 1.78
T2 = 5.13
Tp = 0.5
Gs = tf(K,[T1*T2, T1+T2, 1],'IODelay',T0)
Gz = c2d(Gs,Tp,"zoh")
Kd = dcgain(Gz);
Kd_s = dcgain(Gs);
figure
step(Gs)
hold on
step(Gz)
xlabel("Czas")
ylabel("Wartości na wyjściu y")
title("Odpowiedzi skokowe")
legend("Trasnmitancja ciągła", "Transmitancja dyskretna",Location="southeast")
print('zad1portrans.png','-dpng','-r400')