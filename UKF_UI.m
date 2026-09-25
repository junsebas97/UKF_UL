function [mx, Px, fu] = UKF_UI(y, mx_0, Px_0, fu_0, Q, R, gamma, phi, beta, ...
                               theta)
%{
this function is the UKF-UI filter by Lei et al. It infers the state and
uncertain force of the given system. It uses 4th-order Runge-Kutta for the
time integration and Bouc-Wen model for the restoring force.

INPUTS:
y:     measurements                 [m],  [m/s],   [m/s2]
mx_0:  initial state mean           [m],  [m/s],   [m/s2],  [m]
Px_0:  initial state covariance     [m2], [m2/s2], [m2/s4], [m2]
fu_0:  initial uncertain force      [kN]
Q:     process noise covariance     [m2], [m2/s2], [m2/s4], [m2]
R:     measurement noise covariance [m2], [m2/s2], [m2/s4]
gamma: spread parameter UT          [-]
phi:   spread parameter UT          [-]
beta:  non-Gaussianity parameter UT [-]
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

OUTPUTS:
mx: state mean       [m],  [m/s],   [m/s2],  [m]
Px: state covariance [m2], [m2/s2], [m2/s4], [m2]
fu: uncertain inputs [kN]

BIBLIOGRAPHY: {1} A novel unscented Kalman filter for recursive state-input-
                  system identification of nonlinear systems - Lei Y et al
              {2} A nonlinear Bayesian filter for structural systems with
                  uncertain forces - Delgado Trujillo et al
              {3} Numerical Methods for Engineers (6th edition) - Chapra
                  SC, Canale RP
%}
%% PARAMETRISATION:
% extract the system parameters
f        = theta{1};          % known forces                     [kN]
eta      = theta{2};          % influence vector known loads     [-]
eta_u    = theta{3};          % influence vector uncertain loads [-]
M        = theta{4};          % mass matrix                      [1e3 kg]
C        = theta{5};          % damping matrix                   [kN-s/m]
ke       = theta{6};          % stiffnesses                      [kN/m]
alpha_BW = theta{7};          % post-yield stiffness ratio       [-]
beta_BW  = theta{8};          % parameter of the Bouc-Wen model  [-]
gamma_BW = theta{9};          % parameter of the Bouc-Wen model  [-]
n_BW     = theta{10};         % parameter of the Bouc-Wen model  [-]
dt       = theta{12};         % time step                        [s]
Su       = theta{15};         % displacement selection matrix    [-]
Sv       = theta{16};         % velocity selection matrix        [-]
Sa       = theta{17};         % acceleration selection matrix    [-]

Nt       = size(    y, 2);    % number of times
N_DOFs   = size(    M, 1);    % number of DOFs
N_Su     = size(   Su, 1);    % number of measured displacements
N_Sv     = size(   Sv, 1);    % number of measured velocities
N_Sa     = size(   Sa, 1);    % number of measured accelerations
Nfu      = size(eta_u, 2);    % number of uncertain inputs
Ny       = size(y, 1);        % number of measurements

% remove the acceleration components from the state -- from {2} Page 3
mx_0(2*N_DOFs + (1:N_DOFs))                        = [];   % initial state mean
Px_0(2*N_DOFs + (1:N_DOFs),                     :) = [];   % initial state covar
Px_0(                    :, 2*N_DOFs + (1:N_DOFs)) = [];
Q(   2*N_DOFs + (1:N_DOFs),                     :) = [];   % process noise covar
Q(                       :, 2*N_DOFs + (1:N_DOFs)) = [];

N = size(mx_0, 1);    % number of state components

% compute the weights of the UT
lambda    = phi^2*(N + gamma) - N;                         % [-] -- {1} Page 122
Wm        = NaN(2*N + 1, 1);
Wc        = NaN(1, 2*N + 1);
Wm(1)     = lambda/(N + lambda);                           % [-] -- {2} Eq.6
Wc(1)     = (lambda/(N + lambda)) + (1 - phi^2 + beta);    % [-] -- {2} Eq.7
Wm(2:end) = 1/(2*(N + lambda));                            % [-] -- {2} Eq.8
Wc(2:end) = 1/(2*(N + lambda));                            % [-] -- {2} Eq.8

% define the parameters of the models
theta_g = {M; C; ke; alpha_BW; beta_BW; gamma_BW; n_BW; dt; eta};
theta_h = {M; C; ke; alpha_BW; eta; Su; Sv; Sa};

% calculate influence matrixes of the uncertain inputs
phi_u                               = sparse(N, Nfu);     % state influence
phi_u(N_DOFs + (1:N_DOFs), :)       = M\eta_u;            % [1/1e3kg]
                                                          % -- from {1} Page 123

lambda_u                            = sparse(Ny, Nfu);    % measurements
lambda_u(N_Su + N_Sv + (1:N_Sa), :) = Sa*(M\eta_u);       % influence [1/1e3kg]

%% MAIN:
mx = NaN(  N,    Nt);      mx(:,    1) = mx_0;
fu = NaN(Nfu,    Nt);      fu(:,    1) = fu_0;
Px = NaN(  N, N, Nt);      Px(:, :, 1) = Px_0;

for k = 1:(Nt - 1)
    % for each time step,
    % 1) perform the sigma point calculation step -- {1} Eq.17
    try       sqrt_P = chol((N + lambda)*Px(:, :, k), 'lower');
    catch;    break
    end

    X_kk                         = repmat(mx(:, k), 1, 2*N + 1);    % [m],
    X_kk(:, 1 +           (1:N)) = mx(:, k) + sqrt_P;               % [m/s],
    X_kk(:, 1 + ((N + 1):(2*N))) = mx(:, k) - sqrt_P;               % [m]

    % 2) perform the time prediction step
    X_kp1k       = NaN(N,  2*N + 1);
    y_kp1kp1     = NaN(Ny, 2*N + 1);

    for i = 1:(2*N + 1)
        % for each point
        % 1) evaluate the state
        X_kp1k(:, i) = intgrl_g(X_kk(:, i),         ...      % [m], -- {1} Eq.19
                                f(:, [k, k + 1]),   ...      % [m/s],
                                theta_g)          + ...      % [m]
                       phi_u*fu(:, k)*dt;

        % 2) compute the measurement
        y_kp1kp1(:, i) = h(X_kp1k(:, i), f(:, k + 1), ...    % [m], -- {1} Eq.22
                           theta_h);                         % [m/s],
    end                                                      % [m/s2]

    mx_kp1k = X_kp1k*Wm;                    % predicted    [m],  -- {1} Eq.20
                                            % state vector [m/s],
                                            %              [m]
    X_dist  = X_kp1k - mx_kp1k; 
    Px_kp1k = (Wc.*X_dist)*X_dist' + Q;     % error        [m2], -- {1} Eq.21
                                            % covariance   [m2/s2],
                                            % matrix       [m2]
                                            
    Y_kp1kp1 = y_kp1kp1*Wm;                 % estimated    [m],  -- {1} Eq.22
                                            % measurement  [m/s],
                                            % vector       [m/s2]
    Y_dist   = y_kp1kp1 - Y_kp1kp1;
    Py_kp1   = (Wc.*Y_dist)*Y_dist' + R;    % error        [m2], -- {1} Eq.23
                                            % covariance   [m2/s2],
                                            % matrix       [m2/s4]

    Pxy_kp1  = (Wc.*X_dist)*Y_dist';        % cross-             -- {1} Eq.24
                                            % covariance
                                            % matrix

    % 3) perform the measurement update step
    % 3.1) update the state
    K_kp1        = Pxy_kp1/Py_kp1;                   % Kalman gain  -- {1} Eq.13
    mx(:, k + 1) = mx_kp1k + ...                     % updated [m], -- {1} Eq.27
                   K_kp1*(y(:, k + 1) - Y_kp1kp1);   % state   [m/s],
                                                     % vector  [m]
    % 3.2) estimate the input
    Delta = @(fu) y(:, k + 1)                           - ... % estimation   [m],
                  h(mx(:, k + 1), f(:, k + 1), theta_h) - ... % error        [m/s],
                  lambda_u*fu;                                % -- {1} Eq.29 [m/s2]

    options = optimoptions('lsqnonlin',                        ...
                           'Algorithm', 'levenberg-marquardt', ...
                           'Display',   'off');

    fu(:, k + 1) = lsqnonlin(@(fu) sum(Delta(fu).^2), ...  % uncertan input [kN]
                             fu(:, k), [], [], options);   % -- {1} Page 125

    % 3.3) reevaluate the measurements and reupdate the state
    Y_kp1kp1     = y_kp1kp1*Wm + ...                % estimat [m],  -- {1} Eq.25
                   lambda_u*fu(:, k + 1);           % measure [m/s],
                                                    %         [m/s2]

    mx(:, k + 1) = mx_kp1k + ...                    % updated [m],  -- {1} Eq.27
                   K_kp1*(y(:, k + 1) - Y_kp1kp1);  % state   [m/s],
                                                    % vector  [m]

    Px_aux          = Px_kp1k - ...                 % updated [m2], -- {1} Eq.12
                      K_kp1*Py_kp1*K_kp1';          % error   [m2/s2],
    Px(:, :, k + 1) = (1/2)*(Px_aux + Px_aux');     % covaria [m2]
end

% allocate space for acceleration in the state estimates to match the outside
% format -- Page 3
mx_org = mx;
Px_org = Px;
idx    = [1:2*N_DOFs, 3*N_DOFs + (1:N_DOFs)];

mx              = NaN(4*N_DOFs, Nt);              % mean
mx(idx,      :) = mx_org;                         % [m],  [m/s],   [m/s2],  [m]
Px              = NaN(4*N_DOFs, 4*N_DOFs, Nt);    % covariance
Px(idx, idx, :) = Px_org;                         % [m2], [m2/s2], [m2/s4], [m2]
end

function [xdot] = g(x, f, theta)
%{
this is the state model for the UKF-UI

INPUTS:
x:     system state      [m], [m/s], [m]
f:     force             [kN]
theta: system parameters

OUTPUTS:
xdot: system dynamics [m/s], [m/s2], [m/s]
%}
%% PARAMETRISATION:
% get the system parameters
M      = theta{1};      % mass matrix                  [1e3 kg]
C      = theta{2};      % damping matrix               [kN-s/m]
k      = theta{3};      % stiffnesses                  [kN/m]
alpha  = theta{4};      % post-yield stiffness ratio   [-]
beta   = theta{5};      % BW-model parameter           [-]
gamma  = theta{6};      % BW-model parameter           [-]
n      = theta{7};      % BW-model parameter           [-]
eta    = theta{9};      % influence matrix known force [-]

N_DOFs = size(M, 1);    % number of DOFs               [-]

% extract the state -- {1} Page 122
u       = x(             1:    N_DOFs);    % displacements    [m]
v       = x((  N_DOFs + 1):(2*N_DOFs));    % velocities       [m/s]
xi      = x((2*N_DOFs + 1):       end);    % BW displacements [m]

vi_vim1 = diff([0; v]);                    % velocity drifts  [m/s]

%% MAIN:
% evaluate the restoring force
F = alpha.*k.*diff([0; u]) + (1 - alpha).*k.*xi;    % spring's [kN] -- {2} Eq.50
F = F - [F(2:end); 0];                              % DOFs's   [kN]

% compute the dynamics -- {1} Eqs.2 and 31
xdot = [v;                                                              % [m/s]
        M\(eta*f - C*v - F);                                            % [m/s2]
        vi_vim1 - beta.*abs(vi_vim1).*(abs(xi).^(n - 1)).*xi - ...      % [m/s]
                                      gamma.*vi_vim1.*(abs(xi).^n)];
end

function [x_kp1] = intgrl_g(x_k, f, theta)
%{
this function is the 4th-order Runge-Kutta method

INPUTS:
x:     initial state          [m], [m/s], [m]
f:     initial and end forces [kN]
theta: system parameters

OUTPUTS:
x_kp1: new state [m], [m/s], [m]
%}

dt = theta{8};             % time step    [s]
df = f(:, 2) - f(:, 1);    % input change [kN]

% calculate the slopes -- {3} Eqs.25.40a to 25.40d
k1 = g(x_k,               f(:, 1),            theta);    % [m/s], [m/s2], [m/s]
k2 = g(x_k + (1/2)*k1*dt, f(:, 1) + (1/2)*df, theta);    % [m/s], [m/s2], [m/s]
k3 = g(x_k + (1/2)*k2*dt, f(:, 1) + (1/2)*df, theta);    % [m/s], [m/s2], [m/s]
k4 = g(x_k +       k3*dt, f(:, 1) +       df, theta);    % [m/s], [m/s2], [m/s]

x_kp1 = x_k + ...                           % new state [m], -- {3} Eq.25.40
       (1/6)*(k1 + 2*k2 + 2*k3 + k4)*dt;    %           [m/s]
                                            %           [m]
end

function [y] = h(x, f, theta)
%{
this is the measurement model for the UKF-UI

INPUTS:
x:     system state      [m], [m/s], [m]
f:     force             [kN]
theta: system parameters

OUTPUTS:
y: measurements [m], [m/s], [m/s2]
%}

% get the system parameters and state
M      = theta{1};      % mass matrix                   [1e3 kg]
C      = theta{2};      % damping matrix                [kN-s/m]
k      = theta{3};      % stiffnesses                   [kN/m]
alpha  = theta{4};      % post-yield stiffness ratio    [-]
eta    = theta{5};      % influence matrix known forces [-]
Su     = theta{6};      % displacement selection matrix [-]
Sv     = theta{7};      % displacement selection matrix [-]
Sa     = theta{8};      % displacement selection matrix [-]

N_DOFs = size(M, 1);    % number of DOFs                [-]

% extrat the state -- {1} Page 122
u       = x(             1:    N_DOFs);    % displacements    [m]
v       = x((  N_DOFs + 1):(2*N_DOFs));    % velocities       [m/s]
xi      = x((2*N_DOFs + 1):       end);    % BW displacements [m]

% evaluate the restoring force
F = alpha.*k.*diff([0; u]) + (1 - alpha).*k.*xi;    % spring's [kN] -- {2} Eq.50
F = F - [F(2:end); 0];                              % DOFs'    [kN]

% compute the measuremenets
y = [Su*u;                         % [m],
     Sv*v;                         % [m/s],
     Sa*(M\(eta*f - C*v - F))];    % [m/s2]
end