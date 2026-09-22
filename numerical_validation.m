%{
This code implements the nonlinear Bayesian filter with uncertain forces,
and compares it agains the UKF-UI of Lei Y et al. It estimates the states
and the uncertain forces of the system from noisy measurements

MADE BY:      junsebas97
BIBLIOGRAPHY: A nonlinear Bayesian filter for structural systems with
              uncertain forces - Delgado Trujillo et al
%}
close all; clear all; clc
rng(1234)

%% INPUTS:
% system:
m         = 2*ones(6, 1);              % mass                 [1e3 kg]
c         = 0.5*ones(6, 1);            % damping              [kN-s/m]
k         = [10; 10; 10; 10; 5; 5];    % stiffness            [kN/m]
alpha_BW  = 0.3*ones(6, 1);            % stiffness ratio      [-]
beta_BW   = 80*ones(6, 1);             % BW model parameter   [-]
gamma_BW  = 40*ones(6, 1);             % BW model parameter   [-]
n_BW      = 2*ones(6, 1);              % BW model parameter   [-]
x0        = zeros(24, 1);              % initial condition    [m], [m/s], [m/s2], [m]
rho_infty = 0.8;                       % TIA disipation       [-]
tolerance = 1e-3;                      % iterations tolerance [-]
max_iter  = 20;                        % maximum iterations   [-]

% force:
t    = (0:0.05:60)';                           % times                      [s]
S_fk = zeros(6, 1);                            % known's influence matrix   [-]
S_fu = sparse([4, 6], [1, 2], [1, 1], 6, 2);   % uncertn's influence matrix [-]
f_k  = zeros(1, size(t, 1));                   % known forces               [kN]
f_u  = [-0.1*(2*sin(0.4*t') + sin(0.25*t'));   % uncertain forces           [kN]
        0.25*(1 - cos(0.5*t'))];

% measurements:
noise_ratio = 0.05;                                 % noise to signal ratio  [-]
Su          = sparse(        1,               ...   % displacement selection [-]
                             4,         1, 1, 6);
Sv          = sparse(        1,               ...   % velocity selection     [-]
                             4,         1, 1, 6);
Sa          = sparse([1, 2, 3],               ...   % acceleration selection [-]
                     [1, 2, 6], [1, 1, 1], 3, 6);

% filter:
alphaBW_fltr = 0.4*ones(6, 1);    % post-yield stiffness ratio [-]
nBW_fltr     = 2.5*ones(6, 1);    % BW model parameter         [-]

kappa        = 1;                 % spread factor UT           [-]
alpha        = 1;                 % spread factor UT           [-]
beta         = 0;                 % non-Gaussianity factor UT  [-]

mfu_0 = (5e-1)*ones(2, 1);    % initial load mean   [kN]
Pfu_0 = (2e-1)*eye(2);        % initial load covar  [kN2]
E     = (1e-1)*eye(2);        % random walk covar   [kN2]

mx_0  = (1e-1)*randn(24, 1);  % initial state mean  [m],  [m/s],   [m/s2],  [m]
Px_0  = (1e-2)*eye(24);       % initial state covar [m2], [m2/s2], [m2/s4], [m2]
Q     = (1e-6)*eye(24);       % process noise covar [m2], [m2/s2], [m2/s4], [m2]

R     = diag([1e-2, ...       % covariance          [m2],
              1e-2, ...       % measurement noise   [m2/s2],
              1e-2, ...       %                     [m2/s4],
              1e-2, ...       %                     [m2/s4],
              1e-1]);         %                     [m2/s4]

%% PARAMETERS:
dt     = t(2) - t(1);       % time step [s]
N_DOFs = size(    m, 1);    % number of DOFs
n_fu   = size(mfu_0, 1);    % dimensionality of uncertain loads
nx     = size( mx_0, 1);    % dimensionality of the state

% define the system parameters
M       = sparse(diag(m));           % mass matrix    [1e3 kg]
C       = sparse(N_DOFs, N_DOFs);    % damping matrix [kN-s/m]
C(1, 1) = c(1);

for i = 2:N_DOFs
    C([i - 1, i], [i - 1, i]) = C([i - 1, i], [i - 1, i]) + [ c(i), -c(i);
                                                             -c(i),  c(i)];
end

theta = {f_k; S_fk; S_fu; M; C; k; alphaBW_fltr; beta_BW; gamma_BW; nBW_fltr;
         rho_infty; dt; tolerance; max_iter; Su; Sv; Sa};

%% MAIN:
% create target data (measurements and states)
f_data = S_fk*f_k + S_fu*f_u;    % [kN] -- Eq.13

[x, y] = get_data(x0, f_data, M, C, k, alpha_BW, beta_BW, gamma_BW, n_BW, ...
                  rho_infty, dt, tolerance, max_iter, Su, Sv, Sa, noise_ratio);

% apply the filters
mx  = cell(3, 1);
Px  = cell(3, 1);
mfu = cell(3, 1);
Pfu = cell(3, 1);
[mx{1}, Px{1}, mfu{1}, Pfu{1}] = UKF_UL(y, mx_0, Px_0, mfu_0, Pfu_0, E, Q,  ...
                                        R, kappa, alpha, beta, theta);

[mx{2}, Px{2}, mfu{2}, Pfu{2}] = AUKF(y, mx_0, Px_0, mfu_0, Pfu_0, E, Q, R, ...
                                      kappa, alpha, beta, theta);

[mx{3}, Px{3}, mfu{3}]         = UKF_UI(y, mx_0, Px_0, mfu_0, Q, R, kappa,  ...
                                        alpha, beta, theta);

% evaluate the RMSE of the estimates of the states and the uncertain loads
xRMSE = NaN(  nx, 3);
fRMSE = NaN(n_fu, 3);
for i = 1:3
    xRMSE(:, i) = rmse( mx{i},   x, 2);
    fRMSE(:, i) = rmse(mfu{i}, f_u, 2);
end

%{
NOTE:
The units are as follows
    x, mx, and xRMSE ---> [m],  [m/s],   [m/s2],  [m]
    Px               ---> [m2], [m2/s2], [m2/s4], [m2]
    y                ---> [m],  [m/s],   [m/s2]
    mfu and fRMSE    ---> [kN]
    Pfu              ---> [kN2]
%}

%% REPORT:
print_RMSE(xRMSE, ["UKF-UL", "AUKF", "UKF-UI"])
print_RMSE(fRMSE, ["UKF-UL", "AUKF", "UKF-UI"])

compare_estimation(  x,  y, t,  mx, {Su; Sv; Sa; sparse(0, 6)}, ["UKF-UL", "AUKF", "UKF-UI"])
compare_estimation(f_u, [], t, mfu,             {sparse(0, 2)}, ["UKF-UL", "AUKF", "UKF-UI"])

plot_estimation(  x,  y, t,  mx{1},  Px{1}, {Su; Sv; Sa; sparse(0, 6)})
plot_estimation(f_u, [], t, mfu{1}, Pfu{1},             {sparse(0, 2)})