% Czyszczenie przestrzeni roboczej i okna poleceń
clear
clc

% Warunki początkowe do edycji

delta_t = 0.001; %krok symulacji [s]

m = 1; % masa [kg]
r = 0.2;  % promień [m]
C = 0.25; % Współczynnik oporu

v = 100; % prędkość początkowa [m/s]
alfa = 45;  % kąt rzutu [°]
h = 0.001;  % początkowa wysokość [m]

% Parametry wiatru
wind10 = 5;   % prędkość wiatru na wysokości 10m
wind_h10 = 10;  % wysokość pomiaru - standardowo 10m
terrain_alfa = 0.2; % współczynnik szorstkości terenu
wind_angle = 0;    % kąt wiatru

% koniec bloku do edycji


%Ustawianie początkowych wartości
x(1) = 0;
y(1) = h;
t(1) = 0;

vx(1) = v * cosd(alfa); %wyliczanie składowej poziomej prędkości
vy(1) = v * sind(alfa); %wyliczanie składowej pionowej prędkości 

i = 1; % Ustawienie licznika/indeksu


% Pętla dla ruchu pocisku dopóki nie dotknie ziemi
while (y(i) > 0)

    % Obliczenia gęstości powietrza dla różnych warst atmosfery
    if h < 11000
        L = 0.0065;
        T = 288.5;
    elseif h <= 20000 && h >= 11000
        L = 0;
        T=216.65;
    elseif h>20000
        L = -0.001;
        T = 288.15;
    end
    
    rho = 1.225 * (1 - (L * h) / T)^((0.0341628 / L) - 1); % wyliczanie gęstości


    wind_speed = wind10*(h/wind_h10)^terrain_alfa; % wyliczenie prędkości wiatru na aktualnej wysokości

    % Obliczenia oporu powietrza
    A = (pi * (r * 2)^2) / 4;
    b = (rho * C * A) / (2 * m);
    ax = -(b / m) * vx^2;
    ay = -9.81 - (b / m) * vy^2;

    % Wpływ wiatru
    ax_wind = (wind_speed * cosd(wind_angle));
    ax = ax + ax_wind;

    ay_wind = (wind_speed * cosd(wind_angle));
    ay = ay + ay_wind;

    % Aktualizacja prędkości pocisku
    v = sqrt(vx^2 + vy^2);
    vx = vx + ax * delta_t;
    vy = vy + ay * delta_t;

    % Aktualizacja pozycji i czasu
    x(i + 1) = x(i) + vx * delta_t + 0.5 * ax * delta_t^2;
    y(i + 1) = max(0, y(i) + vy * delta_t + 0.5 * ay * delta_t^2);
    t(i + 1) = t(i) + delta_t;

    % Zapis danych o przyspieszeniu i prędkości do późniejszego wyświetlania
    acceleration_x(i) = ax;
    velocity_x(i) = vx;

    % Zapis informacji o trajektorii
    max_height = max(y);
    total_distance = x(end);
    flight_time = t(end);

    % Aktualizacja wysokości dla kolejnej iteracji
    h = h - vy * delta_t;

    i = i + 1;
end


figure;

% Wykres toru lotu pocisku
subplot(3, 1, 1);
plot(x, y, 'b');
xlabel('Dystans poziomy (m)');
ylabel('Dystans pionowy (m)');
title('Tor lotu pocisku');
axis equal;

% Wykres zmian przyspieszenia poziomego
subplot(3, 1, 3);
plot(t(1:i-1), acceleration_x, 'r');
xlabel('Czas (s)');
ylabel('Ax (m/s^2)');
title('Zmiany w przyspieszeniu poziomym');


% Wykres zmian prędkości poziomej
subplot(3, 1, 2);
plot(t(1:i-1), velocity_x, 'g');
xlabel('Czas (s)');
ylabel('Vx (m/s)');
title('Zmiany w prędkości poziomej');


fprintf('Dystans w osi x: %.2f m\n', total_distance);
fprintf('Maksymalna wysokość: %.2f m\n', max_height);
fprintf('Czas trwania lotu: %.2f s\n', flight_time);
