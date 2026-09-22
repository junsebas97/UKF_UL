%{
This code implements the UKF-UL. It compares the filter performance with
different sensor layouts.

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
noise_ratio = 0.05;                                               % noise2signal
                                                                  % ratio  [-]
Su          = {sparse(1, 4, 1, 1, 6);                             % displacement
               sparse(1, 4, 1, 1, 6);                             % selection
               sparse(1, 4, 1, 1, 6);                             % matrixes [-]
               sparse(1, 4, 1, 1, 6)};
Sv          = {sparse(1, 4, 1, 1, 6);                             % velocity
               sparse(1, 4, 1, 1, 6);                             % selection
               sparse(0, 6);                                      % matrixes [-]
               sparse(0, 6)};
Sa          = {sparse([1, 2, 3], [1, 2, 6], [1, 1, 1], 3, 6);     % acceleration
               sparse([1, 2],    [1,    6], [1,    1], 2, 6);     % selection
               sparse([1, 2],    [1,    6], [1,    1], 2, 6);     % matrixes [-]
               sparse( 1,                6,         1, 1, 6)};

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

R     = {diag([1e-2, 1e-2, 1e-2, 1e-2, 1e-1]);    % covariance measurement
         diag([1e-2, 1e-2,       1e-2, 1e-1]);    % noise [m2], [m2/s2], [m2/s4]
         diag([1e-2, 1e-2,             1e-1]);
         diag([1e-2,                   1e-1])};

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
         rho_infty; dt; tolerance; max_iter; []; []; []};

%% MAIN:
% define the force of the system
f_data = S_fk*f_k + S_fu*f_u;    % [kN] -- Eq.13

% perform the filtering
N_analy = numel(Su);
y       = cell(N_analy, 1);
mx      = cell(N_analy, 1);
mfu     = cell(N_analy, 1);
Px      = cell(N_analy, 1);
Pfu     = cell(N_analy, 1);
xRMSE   = NaN( nx, N_analy);
fRMSE   = NaN(n_fu, N_analy);
for i = 1:N_analy
    % in each analysis,
    % 1) incorporate the selection matrices into the system parameters
    theta(15:17) = {Su{i}, Sv{i}, Sa{i}};

    % 2) create target data (measurements and states)
    rng(1234)
    [x, y{i}] = get_data(x0, f_data, M, C, k, alpha_BW, beta_BW, gamma_BW, ...
                         n_BW, rho_infty, dt, tolerance, max_iter, Su{i},  ...
                         Sv{i}, Sa{i}, noise_ratio);

    % 3) apply the proposed filter
    [mx{i}, Px{i}, mfu{i}, Pfu{i}] = UKF_UL(y{i}, mx_0, Px_0, mfu_0,  Pfu_0, ...
                                            E, Q, R{i}, kappa, alpha, beta,  ...
                                            theta);
    % 4) evaluate the RMSE
    xRMSE(:, i) = rmse( mx{i},   x, 2);
    fRMSE(:, i) = rmse(mfu{i}, f_u, 2);

    % 5) plot the filter's estimation
    plot_estimation(  x, y{i}, t,  mx{i},  Px{i}, {Su{i}; Sv{i}; Sa{i}; sparse(0, 6)})
    plot_estimation(f_u,   [], t, mfu{i}, Pfu{i},                      {sparse(0, 2)})
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
compare_estimation(x, [], t, mx, {sparse(0, 6); sparse(0, 6); sparse(0, 6); sparse(0, 6)}, ...
                   ["1u-1v-3a", "1u-1v-2a", "1u-0v-2a", "1u-0v-1a"])
compare_estimation(f_u, [], t, mfu, {sparse(0, 2)}, ...
                   ["1u-1v-3a", "1u-1v-2a", "1u-0v-2a", "1u-0v-1a"])

print_RMSE(xRMSE, ["1u-1v-3a", "1u-1v-2a", "1u-0v-2a", "1u-0v-1a"])
print_RMSE(fRMSE, ["1u-1v-3a", "1u-1v-2a", "1u-0v-2a", "1u-0v-1a"])