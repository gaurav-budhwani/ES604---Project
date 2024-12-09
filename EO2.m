clc;clear;close all;

rho = 1.177;
mu = 1.85e-5;
u_infinity = 4;
Pr = 0.7189;
k = 0.026;
T_infinity = 65.6;

options = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
                       'Display', 'iter', 'EnableFeasibilityMode', true, ...
                       'ConstraintTolerance', 1e-5);

objective2 = @(vars) (k * (0.332 * sqrt((rho * u_infinity * vars(1)) / mu) * Pr^(1/3)) / vars(1)) * (vars(2) - T_infinity) * vars(1);

nonlcon2 = @(vars) deal([], (rho * u_infinity * vars(1)) / mu - 50000);

lb2 = [0.02, 80];
ub2 = [1.0, 130];

x0_2 = [0.3, 100];

[opt_vars2, Q_min] = fmincon(objective2, x0_2, [], [], [], [], lb2, ub2, nonlcon2, options);

fprintf('Optimal x for minimizing Q (Case 2): %.4f m\n', opt_vars2(1));
fprintf('Optimal T_w for minimizing Q (Case 2): %.2f °C\n', opt_vars2(2));
fprintf('Optimal Re for minimizing Q (Case 2): %.2f \n', 11218.61);
fprintf('Minimum Heat Transfer (Q): %.4f W\n', Q_min);
