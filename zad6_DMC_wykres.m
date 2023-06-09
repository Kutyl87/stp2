% Wykreślenie wykresu liniowego
plot(T0_prop, K_prop_dmc, 'b-'); % Wykres liniowy (kolor niebieski)

% Wypełnienie obszaru poniżej wykresu liniowego
hold on;
x = [T0_prop, fliplr(T0_prop)];
y = [K_prop_dmc, zeros(size(K_prop_dmc))];
fill(x, y, 'b', 'FaceAlpha', 0.3);
hold off;

% Ustawienia osi i etykiet
xlabel('T_{0}/T_{0}^{nom}');
ylabel('K_{0}/K_{0}^{nom}');
title('Obszary stabilności i niestabilności dla regulacji DMC');
grid on;
legend('Obszar niestabilności', 'Obszar stabilności');