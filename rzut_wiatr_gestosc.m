% Clear workspace and command window
clear
clc

% Constants for Bullet
m = 5;       % kg (mass)
r = 0.2;     % m (radius)
C = 0.25;    % Drag Coefficient of a Sphere

% Initial Conditions
delta_t = 0.001; % s
h = 5;         % m (initial height) czemu tutaj nie może być 0???

% Case 1: Free-Falling Bullet
x(1) = 0;       % m
y(1) = h;       % m (starting from the initial height)
v = 500;        % m/s (velocity)

% Case 2: Projectile Bullet With Drag and Wind
theta2 = 45;    % degree
vx(1) = v * cosd(theta2);
vy(1) = v * sind(theta2);
t(1) = 0;
i = 1;          % Sets Counter/Index

% Wind parameters
wind_speed = 10;    % m/s
wind_angle = 0;    % degree

% Start Loop For Projectile Motion with Air Resistance and Wind
while (y(i) > 0)
    % Calculate gęstość powietrza using the formula mentioned earlier
    if h < 11000
        L = 0.0065; % Example value for variable L, you can replace it with your specific variation
        T = 288.5;
    elseif h <= 20000 && h >= 11000
        L = 0;
        T=216.65;
    elseif h>20000
        L = -0.001;
        T = 288.15;
    end
    
    rho = 1.225 * (1 - (L * h) / T)^((0.0341628 / L) - 1);

    % Drag calculation
    A = (pi * (r * 2)^2) / 4;
    b = (rho * C * A) / (2 * m);
    ax = -(b / m) * vx^2;

    % Wind effect
    ax_wind = (wind_speed * cosd(wind_angle));
    ax = ax + ax_wind;

    % Vertical motion
    ay = -9.81 - (b / m) * vy^2;
    v = sqrt(vx^2 + vy^2);
    vx = vx + ax * delta_t;
    vy = vy + ay * delta_t;

    % Update position and time
    x(i + 1) = x(i) + vx * delta_t + 0.5 * ax * delta_t^2;
    y(i + 1) = max(0, y(i) + vy * delta_t + 0.5 * ay * delta_t^2); % Ensure y doesn't go below 0
    t(i + 1) = t(i) + delta_t;

    % Store acceleration and velocity data for plotting
    acceleration(i) = ax;
    velocity_x(i) = vx;
    velocity_y(i) = vy;
    velocity(i)=v;
    % Update height for the next iteration
    h = h - vy * delta_t;

    i = i + 1;
end

% Plot Projectile Motion with Drag and Wind
subplot(3, 1, 1);
plot(x, y, 'b');
xlabel('Horizontal Distance (m)');
ylabel('Vertical Distance (m)');
title('Projectile Path with Wind');
legend('Case 2 with Wind');

% Plot Changes in Horizontal Acceleration
subplot(3, 1, 2);
plot(t(1:i-1), acceleration, 'r');
xlabel('Time (s)');
ylabel('Ax (m/s^2)');
title('Changes in Horizontal Acceleration');

% Plot Changes in Horizontal Velocity
subplot(3, 1, 3);
plot(t(1:i-1), velocity, 'g');
xlabel('Time (s)');
ylabel('Vx (m/s)');
title('Changes in Horizontal Velocity');

figure;
plot(x(1:i-1), y(1:i-1), 'k');
xlabel('Horizontal Distance (m)');
ylabel('Vertical Distance (m)');
title('Vertical Position vs Horizontal Position');
axis equal;

