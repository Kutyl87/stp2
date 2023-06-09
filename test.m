x = [1 2 3 4 5];  % Przykładowy wektor x
y = [2 4 1 3 5];  % Przykładowy wektor y
stairs(x, y);  % Tworzenie wykresu typu stairs

% Znalezienie indeksu dla punktu x = 3
index = find(x == 2);

% Wysokość schodka w punkcie x = 3
difference = y(index)