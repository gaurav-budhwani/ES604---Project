clc;clear;close all

T_inf = 65.6;
k_air = 0.026;
mu_air = 1.85e-5;
Pr = 0.7189;
rho_air = 1.177;
u_inf = 0.8;
w = 1;

x0 = [0.2, 100];  

objective = @(vars) 0.977 * 4.64 * vars(1) / sqrt((rho_air * u_inf * vars(1)) / mu_air) * Pr^(1/3);

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

[opt_vars, opt_delta_t] = fmincon(objective, x0, [], [], [], [], [0.2, 80], [1, 130], @constraints, options);

opt_x = opt_vars(1);
opt_T = opt_vars(2);
opt_Re = (rho_air * u_inf * opt_x) / mu_air;
opt_Nu = 0.332 * sqrt(opt_Re) * Pr^(1/3);
opt_h = opt_Nu * k_air / opt_x;
opt_Q = 2 * opt_h * opt_x * w * (opt_T - T_inf);

fprintf('Optimal x: %.4f\n', opt_x);
fprintf('Optimal T: %.4f\n', opt_T);
fprintf('Optimal Re: %.4f\n', opt_Re);
fprintf('Optimal Q: %.4f\n', opt_Q);
fprintf('Optimal delta_t: %.4f\n', opt_delta_t);

function [c, ceq] = constraints(vars)
    x = vars(1);
    T = vars(2);
    T_inf = 65.6;
    k_air = 0.026;
    mu_air = 1.85e-5;
    Pr = 0.7189;
    rho_air = 1.177;
    u_inf = 0.8;
    w = 1;
    
    Re = (rho_air * u_inf * x) / mu_air;
    Nu = 0.332 * sqrt(Re) * Pr^(1/3);
    h = Nu * k_air / x;
    Q = 2 * h * x * w * (T - T_inf);
    
    c = [
        -Re;                 % Re >= 1
        Re - 50000;          % Re <= 50000
        140 - Q;             % Q >= 140
        Q - 190;             % Q <= 190
        0.2 - x;             % x >= 0.2
        x - 1;               % x <= 1
        80 - T;              % T >= 80
        T - 130;             % T <= 130
    ];
    
    ceq = [];
end
