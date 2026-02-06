clear; close all; 
%Just a toy demo for the point particle MSD

%% Physical parameters
kB = 1;       % Boltzmann constant (arbitrary units)
T = 1;        % Temperature
eta = 1;      % Fluid viscosity
a = 1;        % Sphere radius

%% Simulation parameters
D = kB * T / (6 * pi * eta * a);   % Theoretical diffusion coefficient
dt = 1e-4;                         % Time step
Nsteps = 1e5;                      % Total number of steps
d = 2;                             % Dimension

%% Preallocate
X = zeros(Nsteps, d);             % Store positions
x = zeros(1, d);                  % Initial position

%% Simulate Brownian motion using Euler-Maruyama
for i = 2:Nsteps
    dx = sqrt(2 * D * dt) * randn(1, d);
    x = x + dx;
    X(i, :) = x;
end

%% Compute MSD
lags = round(logspace(1, 5, 100));  % Logarithmically spaced lags
max_lag = floor(Nsteps/10);
lags = unique(lags(lags <= max_lag));
msd = zeros(size(lags));
for j = 1:length(lags)
    tau = lags(j);
    displacements = X(1+tau:end,:) - X(1:end-tau,:);
    sq_disp = sum(displacements.^2, 2);
    msd(j) = mean(sq_disp);
end

%% Plot MSD
figure;
loglog(dt * lags, msd, 'b', 'LineWidth', 2);
hold on;
loglog(dt * lags, 2 * d * D * dt * lags, 'k--', 'LineWidth', 2);
xlabel('Lag time \tau'); ylabel('MSD(\tau)');
legend('Simulated MSD', 'Expected 2dD\tau');
title('Brownian motion of a sphere: MSD test');
grid on;

% Choose a range of lags where MSD grows linearly
last_lags = 40;
range = lags(end-last_lags:end);        % last few lags (adjust if needed)
tau_range = dt * range;
msd_range = msd(end-last_lags:end);

% Linear fit
coeffs = polyfit(tau_range, msd_range, 1);
slope = coeffs(1);
D_measured = slope / (2 * d);

fprintf("Measured D = %.5f\n", D_measured);
fprintf("Theoretical D = %.5f\n", D);
fprintf("Relative error = %.2f%%\n", 100 * abs(D_measured - D) / D);

figure()
tau = dt * lags;
plot(tau, msd ./ (2 *d * tau), 'r');
hold on
plot(tau,ones(size(tau))*1/6/pi,'k--')
xlabel('\tau'); ylabel('MSD(\tau) / 6\tau');
title('Extracted D from MSD slope should be 1/6pi');
