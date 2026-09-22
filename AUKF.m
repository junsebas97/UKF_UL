function [mx, Px, mfu, Pfu] = AUKF(y, mx_0, Px_0, mfu_0, Pfu_0, E, Q, R, ...
                                   kappa, alpha, beta, theta)
%{
this function is the AUKF. It augments the state with uncertain loads and
infers the augmented state with the UKF.

INPUTS:
y:     measurements                             [m],  [m/s],   [m/s2]
mx_0:  initial state mean                       [m],  [m/s],   [m/s2],  [m]
Px_0:  initial state covariance                 [m2], [m2/s2], [m2/s4], [m2]
mfu_0: mean of the initial uncertain load       [kN]
Pfu_0: covariance of the initial uncertain load [kN2]
E:     covariance of uncertain load random walk [kN2]
Q:     covariance process noise                 [m2], [m2/s2], [m2/s4], [m2]
R:     covariance measurement noise             [m2], [m2/s2], [m2/s4]
kappa: spread parameter UT                      [-]
alpha: spread parameter UT                      [-]
beta:  non-Gaussianity parameter UT             [-]
theta: system parameters. They must organized as
    ---                                       ---
   | known forces                      [kN]     |
   | influence vector known loads      [-]      |
   | influence vector uncertain loads  [-]      |
   | mass matrix                       [1e3 kg] |
   | damping matrix                    [kN-s/m] |
   | stiffnesses                       [kN/m]   |
   | post-yield stiffness ratio        [-]      |
   | parameter of the Bouc-Wen model   [-]      |
   | parameter of the Bouc-Wen model   [-]      |
   | parameter of the Bouc-Wen model   [-]      |
   | high energy dissipation           [-]      |
   | time step                         [s]      |
   | tolerance Newton-Raphson          [-]      |
   | maximum iterations Newton-Raphson [-]      |
   | displacement selection matrix     [-]      |
   | velocity selection matrix         [-]      |
   | acceleration selection matrix     [-]      |
   ---                                        ---

OUTPUTs:
mx:  state mean                    [m],  [m/s],   [m/s2],  [m]
Px:  state covariance              [m2], [m2/s2], [m2/s4], [m2]
mfu: mean of uncertain loads       [kN]
Pfu: covariance of uncertain loads [kN2]

MADE BY:      junsebas97
BIBLIOGRAPHY: {1} A nonlinear Bayesian filter for structural systems with
                  uncertain forces - Delgado Trujillo et al
              {2} A Time Integration Algorithm for Structural Dynamics With
                  Improved Numerical Dissipation: The Generalized-alpha
                  Method - Chung J, Hulbert G
%}
%% PARAMETRISATION:
% extract the system parameters
f_k       = theta{1};              % known loads                         [kN]
rho_infty = theta{11};             % high energy dissipation of the TIA  [-]
Su        = theta{15};             % displacement selection matrix       [-]
Sv        = theta{16};             % velocity selection matrix           [-]
Sa        = theta{17};             % acceleration selection matrix       [-]
theta     = theta(1:14);           % system parameters for process model

N_DOFs    = size(   Su, 2);        % number of DOFs
nfu       = size(mfu_0, 1);        % dimensionality of the uncertain loads
N_Su      = size(   Su, 1);        % number displacement measurements
N_Sv      = size(   Sv, 1);        % number velocities measurements
N_Sa      = size(   Sa, 1);        % number acceleration measurements
nt        = size(  f_k, 2);        % number of time steps
nx        = size( mx_0, 1);        % dimensionality of the state
Ny        = N_Su + N_Sv + N_Sa;    % total number of measurements
nz        = nx + nfu;              % augmented dimensionality

% get the measurement matrix
H                                                = zeros(Ny, nx);    % {1} Eq.11
H(              (1:N_Su),            (1:N_DOFs)) = Su;               % [-]
H(       N_Su + (1:N_Sv),   N_DOFs + (1:N_DOFs)) = Sv;               % [-]
H(N_Su + N_Sv + (1:N_Sa), 2*N_DOFs + (1:N_DOFs)) = Sa;               % [-]

% evaluate the known loads at the sampling time
alpha_f    = rho_infty/(rho_infty + 1);          % force factor [-] -- {2} Eq.25
f_k(2:end) = alpha_f*f_k(:, (2:end) - 1) + ...   % known loads  [kN]
             (1 - alpha_f)*f_k(:, (2:end));

% compute the weights of the UT
lambda    = alpha^2*(nz + kappa) - nz;                       % [-] -- {1} Page 3
Wm        = NaN(2*nz + 1, 1);
Wc        = NaN(1, 2*nz + 1);
Wm(1)     = lambda/(nz + lambda);                            % [-] -- {1} Eq.6
Wc(1)     = (lambda/(nz + lambda)) + (1 - alpha^2 + beta);   % [-] -- {1} Eq.7
Wm(2:end) = 1/(2*(nz + lambda));                             % [-] -- {1} Eq.8
Wc(2:end) = 1/(2*(nz + lambda));                             % [-] -- {1} Eq.8

%% MAIN:
% augment the system state with the uncertain inputs
mz = NaN(nz,     nt);
Pz = NaN(nz, nz, nt);

mz(:,    1)        = [mx_0; mfu_0];           % [m],  [m/s],   [m/s2],  [m],  [kN]
Pz(:, :, 1)        = blkdiag(Px_0, Pfu_0);    % [m2], [m2/s2], [m2/s4], [m2], [kN2]
Qz                 = blkdiag(Q, E);           % [m2], [m2/s2], [m2/s4], [m2], [kN2]
H(:, nx + (1:nfu)) = 0;

% apply the UKF to the augmented state
for i = 2:nt
    % for each time step:
    % 1) incorporate the current known inputs into the system parameters
    theta{1} = f_k(:, i);    % [kN]

    % 2) perform the prediction
    % 2.1) create the sigma points of the UT -- {1} Eqs.1 to 3
    sqrt_Pz               = sqrt(nz + lambda)*chol(Pz(:, :, i - 1), 'lower');
    Z                     = NaN(nz, 2*nz + 1);
    Z(:,               1) = mz(:, i - 1);              % [m],
    Z(:, 1 +      (1:nz)) = mz(:, i - 1) + sqrt_Pz;    % [m/s],
    Z(:, 1 + nz + (1:nz)) = mz(:, i - 1) - sqrt_Pz;    % [m/s2],
                                                       % [m],
                                                       % [kN]
    % 2.2) propagate the sigma points
    g_Z = NaN(nz, 2*nz + 1);
    for j = 1:(2*nz + 1)
        g_Z(:, j) = g(Z(:, j), theta);     % [m], [m/s], [m/s2], [m], [kN]
    end

    % 2.3) compute the statistics
    mz_iim1 = g_Z*Wm;                      % [m], [m/s], [m/s2], [m], [kN]
                                           % -- {1} Eq.4
    g_dist   = g_Z - mz_iim1;
    Pz_iim1 = (Wc.*g_dist)*g_dist' + Qz;   % [m2], [m2/s2], [m2/s4], [m2], [kN2]
                                           % -- {1} Eq.5
    % 3) update the augmented state
    my          = H*mz_iim1;                     % [m],  [m/s],   [m/s2]
    Py          = H*Pz_iim1*H' + R;              % [m2], [m2/s2], [m2/s4]

    K           = (Pz_iim1*H')/Py;
    mz(:,    i) = mz_iim1 + K*(y(:, i) - my);    % [m],  [m/s],   [m/s2],  [m],  [kN]
    Pz_aux      = Pz_iim1 - K*Py*K';             % [m2], [m2/s2], [m2/s4], [m2], [kN2]
    Pz(:, :, i) = (1/2)*(Pz_aux + Pz_aux');
end

% separate the augmented state
mx  = mz(1:nx,       :);    % state mean       [m2], [m2/s2], [m2/s4], [m2]
Px  = mz(1:nx, 1:nx, :);    % state covariance [m2], [m2/s2], [m2/s4], [m2]

mfu = mz(nx + (1:nfu),               :);    % uncertain load mean       [kN]
Pfu = mz(nx + (1:nfu), nx + (1:nfu), :);    % uncertain load covariance [kN2]
end

function [z_i] = g(z_im1, theta)
%{
this function is the SSTIA-built process model. It uses the generalised-alpha
method to compute the current augmented state with the previous augmented
state. The system employs the Bouc-Wen model for the restoring force.

INPUTS:
z_im1: previous augmented state [m], [m/s], [m/s2], [m], [kN]
theta: system parameters

OUTPUTs:
z_i: current augmented state [m], [m/s], [m/s2], [m], [kN]

MADE BY:      junsebas97
BIBLIOGRAPHY: {1} A nonlinear Bayesian filter for structural systems with
                  uncertain forces - Delgado Trujillo et al
              {2} A Time Integration Algorithm for Structural Dynamics With
                  Improved Numerical Dissipation: The Generalized-alpha
                  Method - Chung J, Hulbert G
%}
%% PARAMETRISATION:
% get the system properties
fk_i      = theta{1};         % current known loads                  [kN]
S_fk      = theta{2};         % influence vector for known loads     [-]
S_fu      = theta{3};         % influence vector for uncertain loads [-]
M         = theta{4};         % mass matrix                          [1e3 kg]
C         = theta{5};         % damping matrix                       [kN-s/m]
k         = theta{6};         % stiffnesses                          [kN/m]
alpha_BW  = theta{7};         % post-yield stiffness ratio           [-]
beta_BW   = theta{8};         % parameter of the Bouc-Wen model      [-]
gamma_BW  = theta{9};         % parameter of the Bouc-Wen model      [-]
n_BW      = theta{10};        % parameter of the Bouc-Wen model      [-]
rho_infty = theta{11};        % high energy dissipation              [-]
dt        = theta{12};        % time step                            [s]
tolerance = theta{13};        % tolerance Newton-Raphson             [-]
max_iter  = theta{14};        % maximum iterations Newton-Raphson    [-]

nfu       = size(S_fu, 2);    % number of uncertain loads
N_DOFs    = size(M, 1);       % numbers of DOFs

% calculate the parameters of the generalised-alpha method
alpha_f = rho_infty/(rho_infty + 1);            % [-] -- {2} Eq.25
alpha_m = (2*rho_infty - 1)/(rho_infty + 1);    % [-] -- {2} Eq.25
gamma   = 1/2 - alpha_m + alpha_f;              % [-] -- {2} Eq.17
beta    = (1/4)*(1 - alpha_m + alpha_f)^2;      % [-] -- {2} Eq.20

% extract the previous state
u_tim1  = z_im1(           (1:N_DOFs));    % displacement     [m]
v_tim1  = z_im1(  N_DOFs + (1:N_DOFs));    % velocity         [m/s]
a_tim1  = z_im1(2*N_DOFs + (1:N_DOFs));    % acceleration     [m/s2]
xi_tim1 = z_im1(3*N_DOFs + (1:N_DOFs));    % BW displacements [m]
fu_im1  = z_im1(4*N_DOFs + (1:nfu));       % uncertain loads  [kN]

% assess the previous restoring force
r_tim1 = alpha_BW.*k.*diff([0; u_tim1]) + ...    % springs force [kN] 
         (1 - alpha_BW).*k.*xi_tim1;             % -- {1} Eq.50
r_tim1 = r_tim1 - [r_tim1(2:end); 0];            % DOFs force    [kN]

%% MAIN:
% evaluate the uncertain load and the external force vector
fu_i = fu_im1;                   % uncertain load        [kN]
f_i  = S_fk*fk_i + S_fu*fu_i;    % external force vector [kN] -- {1} Eq.12

% initalize the displacements and velocities
u_ti = u_tim1;    % [m]
v_ti = v_tim1;    % [m/s]

% compute the equivalent stiffness matrix and force vector
K_hat = M*((1 - alpha_m)/(beta*dt^2))     +                      ... % [kN/m]
        C*((1 - alpha_f)*gamma/(beta*dt));

f_hat = f_i                                                    - ... % [kN]
        (  M*(1 - 1/(2*beta) + alpha_m/(2*beta))               + ...
           C*(dt*(1 - alpha_f)*(1 - gamma/(2*beta)))  )*a_tim1 - ...
        (  C*(1 - gamma/beta + alpha_f*gamma/beta)             - ...
           M*((1 - alpha_m)/(beta*dt))                )*v_tim1 - ...
        (- M*((1 - alpha_m)/(beta*dt^2))                       - ...
           C*((1 - alpha_f)*gamma/(beta*dt))          )*u_tim1 - ...
        alpha_f*r_tim1;

% find the new nodal displacements with Newton-Raphson iteration
for iter = 1:max_iter
    % in each iteration,
    % 1) calculate the BW response
    [xi_ti, r_ti, dr_du] = BW_model(xi_tim1, v_tim1, v_ti, u_ti, ...   % [m],
                                    k, alpha_BW, beta_BW,        ...   % [kN],
                                    gamma_BW, n_BW, dt);               % [kN/m],

    % 2) assess the tangent stiffness and compute the residual force
    K_t   = K_hat + (1 - alpha_f)*dr_du;                  % [kN/m]
    R_for = f_hat - (K_hat*u_ti + (1 - alpha_f)*r_ti);    % [kN]

    % 3) if tolerance is not met, update the current displacements and
    %    velocities
    if norm(R_for) > tolerance
        duN  = K_t\R_for;                                 % increment    [m]
        u_ti = u_ti + duN;                                % displacement [m]
        v_ti = (gamma/(beta*dt))*(u_ti - u_tim1) + ...    % velocity     [m/s]
               (1 - gamma/beta)*v_tim1           + ... 
               dt*(1 - gamma/(2*beta))*a_tim1;
    else
        break
    end
end

% compute the current accelerations and form the current state vector
a_ti = (1/(beta*dt^2))*(u_ti - u_tim1) - ...    % [m/s2] 
       (1/(beta*dt))*v_tim1            - ...
       (1/(2*beta) - 1)*a_tim1;

z_i  = [u_ti; v_ti; a_ti; xi_ti; fu_i];         % [m], [m/s], [m/s2], [m], [kN]
end