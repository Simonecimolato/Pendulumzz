%% SIMONE CIMOLATO
% June 6th, 2024
% N pendulums with masses at each pendulum's end

clear
close all
clc

%% data
N = 3;
g = 9.81;
dampingCoeffs = 0.3 * ones(N, 1);
tSpan = [0 50];

l = ones(N, 1);
m = ones(N, 1);
Y0 = [pi/2; zeros(N-1, 1); zeros(N, 1)];

%% ode
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-7);
[t, y] = ode113(@(t, y) motion(y, g, dampingCoeffs, l, m, N), tSpan, Y0, options);
y(:, 1:N) = wrapToPi(y(:, 1:N));

%% plots
figure
grid on

% position vs time
subplot(2, 2, 1);
p1 = gobjects(N, 1);
hold on

for i = 1:N
    p1(i) = plot(t(1), rad2deg(y(1, i)), LineWidth = 1.5, DisplayName = sprintf('$\\theta_%d$', i));
end

grid on
xlim([min(t) max(t)]);
ylim([min(rad2deg(y(:, 1:N)), [], "all") max(rad2deg(y(:, 1:N)), [], "all")]);
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Angular positions $\theta_i [^\circ]$', 'Interpreter', 'latex')
title('$\theta_i$ vs Time', 'Interpreter', 'latex')

% quiver animation
subplot(2, 2, 2);
p2 = gobjects(N, 1);

xi = zeros(N, 1);
yi = zeros(N, 1);
ui = zeros(N, 1);
vi = zeros(N, 1);

% Components of the first quiver arrow
ui(1) = l(1) * sin(y(1, 1));
vi(1) = -l(1) * cos(y(1, 1));

p2(1) = quiver(xi(1), yi(1), ui(1), vi(1), LineWidth = 1.5, AutoScale = "off");
hold on

[xi, yi, ui, vi] = calcCoord(y, N, l, xi, yi, ui, vi, 1);

for i = 2:N
    p2(i) = quiver(xi(i), yi(i), ui(i), vi(i), LineWidth = 1.5, AutoScale = "off");
end

grid on
axis equal
xlim([-sum(l) sum(l)]);
ylim([-sum(l) sum(l)]);
xlabel('x [m]', 'Interpreter', 'latex')
ylabel('y [m]', 'Interpreter', 'latex')
title('Pendulum animation', 'Interpreter', 'latex')

% speed vs time
subplot(2, 2, 3);
p3 = gobjects(N, 1);
hold on

for i = N+1:2*N
    p3(i - N) = plot(t(1), rad2deg(y(1, i)), LineWidth = 1.5, DisplayName = sprintf('$$\\dot{\\theta_%d}$$', i - N));
end

grid on
xlim([min(t) max(t)]);
ylim([min(rad2deg(y(:, N+1:2*N)), [], "all") max(rad2deg(y(:, N+1:2*N)), [], "all")]);
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Angular speed $\dot{\theta_i} [\frac{^\circ}{s}]$', 'Interpreter', 'latex')
title('$\dot{\theta_i}$ vs Time', 'Interpreter', 'latex')

% phase space
subplot(2, 2, 4);
p4 = gobjects(N, 1);
hold on

for i = 1:N
    p4(i) = plot(rad2deg(y(1, i)), rad2deg(y(1, N+i)), LineWidth=1.5, DisplayName = sprintf('$$\\theta_%d, \\dot{\\theta_%d}$$', i, i));
end

grid on
xlim([min(rad2deg(y(:, 1:N)), [], "all") max(rad2deg(y(:, 1:N)), [], "all")]);
ylim([min(rad2deg(y(:, N+1:2*N)), [], "all") max(rad2deg(y(:, N+1:2*N)), [], "all")]);
xlabel('$\theta_i [^\circ]$', 'Interpreter', 'latex')
ylabel('$\omega_i [\frac{^\circ}{s}]$', 'Interpreter', 'latex')
title('Phase space', 'Interpreter', 'latex')

% Animation loop
for i = 1:length(t)
    subplot(2, 2, 1);
    for j = 1:N
        set(p1(j), 'XData', t(1:i), 'YData', rad2deg(y(1:i, j)));
    end

    subplot(2, 2, 2);
    [xi, yi, ui, vi] = calcCoord(y, N, l, xi, yi, ui, vi, i);
    for j = 1:N
        set(p2(j), 'XData', xi(j), 'YData', yi(j), 'UData', ui(j), 'VData', vi(j));
    end
     
    subplot(2, 2, 3)
    for j = 1:N
        set(p3(j), 'XData', t(1:i), 'YData', rad2deg(y(1:i, j + N)));
    end

    subplot(2, 2, 4)
    for j = 1:N
        set(p4(j), 'XData', rad2deg(y(1:i, j)), 'YData', rad2deg(y(1:i, j + N)));
    end   

    drawnow
end

subplot(2, 2, 1)
legend show
legend('Interpreter','latex')

subplot(2, 2, 3)
legend show
legend('Interpreter','latex')

subplot(2, 2, 4)
legend show
legend('Interpreter','latex')

%% functions definitions
function dydt = motion(y, g, dampCoeffs, lengths, masses, N)
    thetas = y(1:N);
    theta_dots = y(N+1:2*N);
    
    M = zeros(N, N);
    C = zeros(N, 1);
    
    % Compute the mass matrix M
    for i = 1:N
        for j = 1:N
            if i == j
                for k = i:N
                    M(i, j) = M(i, j) + masses(k) * lengths(i)^2;
                end
            else
                for k = max(i, j):N
                    M(i, j) = M(i, j) + masses(k) * lengths(i) * lengths(j) * cos(thetas(i) - thetas(j));
                end
            end
        end
    end
    
    for i = 1:N
        % Gravity term
        for k = i:N
            C(i) = C(i) - masses(k) * g * lengths(i) * sin(thetas(i));
        end
        
        % Damping term
        C(i) = C(i) - dampCoeffs(i) * theta_dots(i);
        
        % Coriolis term
        for j = 1:N
            if i ~= j
                for k = max(i, j):N
                    C(i) = C(i) - masses(k) * lengths(i) * lengths(j) * theta_dots(j)^2 * sin(thetas(i) - thetas(j));
                end
            end
        end
    end
    
    % Solve for theta double dot (angular accelerations)
    theta_ddot = M \ C;
    
    dydt = zeros(2 * N, 1);
    dydt(1:N) = theta_dots;
    dydt(N+1:2*N) = theta_ddot;
end

function [x, y, u, v] = calcCoord(data, N, l, x, y, u, v, tInd)
% Computes the coordinates and quiver components
    if N == 1
        angle = data(tInd, 1);

        u = l * sin(angle);
        v = -l * cos(angle);
    else
        for i = 1:N-1
            angle = sum(data(tInd, 1:i));
                
            % Components of the quiver arrows
            u(i) = l(i) * sin(angle);
            v(i) = -l(i) * cos(angle);
    
            % Base of the subsequent pendulums
            x(i+1) = x(i) + u(i);
            y(i+1) = y(i) + v(i);
        end
        angle = angle + data(tInd, i+1);

        u(i+1) = l(i+1) * sin(angle);
        v(i+1) = -l(i+1) * cos(angle);
    end
end
