
clc;
clear all;
close all;

% Given parameter values
a = 1;
b = 3;
c = 1;
d = 5;
I = 1.8;

gc = 16;
ge = 1.5;
theta = -0.25;
tc = 10;
te = 2;
Vs = 2;
lambda = 10;

r = 0.0068;
xr = 1.6;
s = 4;

% Define the equations
myEquations = @(t, vars) [
    vars(2) - a * vars(1)^3 + b * vars(1)^2 + I - vars(3) + ge * (vars(1) - vars(1, 3)) + (gc * (vars(1) - Vs) / (1 + exp(-lambda * vars(1) * (t - tc) - theta)));
    c - d * vars(1)^2 - vars(2);
    r * (s * (vars(1) - xr) - vars(3))
];

% Set the initial conditions and time span for the integration
initialConditions = [0; 0; 0];
timeSpan = [0 10];

% Solve the system of equations
[t, vars] = ode45(myEquations, timeSpan, initialConditions);

% Plot the variables
figure;
plot(t, vars(:, 1), 'b', 'LineWidth', 2); hold on;
plot(t, vars(:, 2), 'r', 'LineWidth', 2);
plot(t, vars(:, 3), 'g', 'LineWidth', 2);
xlabel('Time');
ylabel('Variable Value');
legend('x', 'y', 'z');
title('Variable Evolution');

% Find equilibrium points
equilibriumPoints = fsolve(@(vars) myEquations(0, vars), initialConditions);

disp('Equilibrium Points:');
disp(equilibriumPoints);

% Calculate the Jacobian matrix at equilibrium points
JacobianMatrix = jacobian(sym(myEquations), sym('vars'));
JacobianAtEquilibrium = subs(JacobianMatrix, sym('vars'), equilibriumPoints);

disp('Jacobian Matrix at Equilibrium:');
disp(JacobianAtEquilibrium);
