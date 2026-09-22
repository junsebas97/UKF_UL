function [x, y] = get_data(x0, f, M, C, k, alpha_BW, beta_BW, gamma_BW,      ...
                           n_BW, rho_infty, dt, tolerance, max_iter, Su, Sv, ...
                           Sa, noise_ratio)
%{
this function creates the data. It returns the noisy measurements, and the
true states. It uses the generalised-alpha method for the time integration
and the Bouc-Wen model for the restoring force.

INPUTS:
x0:          initial system state                      [m], [m/s], [m/s2], [m]
f:           external forces                           [kN]
M:           mass matrix                               [1e3 kg]
C:           damping matrix                            [kN-s/m]
k:           stiffnesses                               [kN/m]
alpha_BW:    post-yield stiffness ratio                [-]
beta_BW:     parameter of the Bouc-Wen model           [-]
gamma_BW:    parameter of the Bouc-Wen model           [-]
n_BW:        parameter of the Bouc-Wen model           [-]
rho_infty:   high energy dissipation                   [-]
dt:          time step                                 [s]
tolerance:   tolerance Newton-Raphson                  [-]
max_iter:    maximum iterations Newton-Raphson         [-]
Su:          displacement selection matrix             [-]
Sv:          velocity selection matrix                 [-]
Sa:          acceleration selection matrix             [-]
noise_ratio: noise-to-signal ratio of the measurements [-]

OUTPUTS:
x: system states [m], [m/s], [m/s2], [m]
y: measurements  [m], [m/s], [m/s2]

MADE BY:      junsebas97
BIBLIOGRAPHY: {1} A nonlinear Bayesian filter for structural systems with
                  uncertain forces - Delgado Trujillo et al
              {2} A Time Integration Algorithm for Structural Dynamics With
                  Improved Numerical Dissipation: The Generalized-alpha
                  Method - Chung J, Hulbert G
%}

%% PARAMETRISATION:
N_DOFs = size( M, 1);           % number of DOFs
N_Su   = size(Su, 1);           % number displacement measurements
N_Sv   = size(Sv, 1);           % number velocities measurements
N_Sa   = size(Sa, 1);           % number acceleration measurements
Nt     = size( f, 2);           % number of times
nx     = 4*N_DOFs;              % dimensionality of the state
Ny     = N_Su + N_Sv + N_Sa;    % total number of measurements

% calculate the parameters of the TIA
alpha_f = rho_infty/(rho_infty + 1);            % [-] -- {2} Eq.25
alpha_m = (2*rho_infty - 1)/(rho_infty + 1);    % [-] -- {2} Eq.25
gamma   = 1/2 - alpha_m + alpha_f;              % [-] -- {2} Eq.17
beta    = (1/4)*(1 - alpha_m + alpha_f)^2;      % [-] -- {2} Eq.20

% get the measurement matrix
H                                                = zeros(Ny, nx);    % {1} Eq.11
H(              (1:N_Su),            (1:N_DOFs)) = Su;               % [-]
H(       N_Su + (1:N_Sv),   N_DOFs + (1:N_DOFs)) = Sv;               % [-]
H(N_Su + N_Sv + (1:N_Sa), 2*N_DOFs + (1:N_DOFs)) = Sa;               % [-]

% extract the previous state -- {1} Page 3
u0  = x0(             1:  N_DOFs);    % displacement    [m]
v0  = x0(  (N_DOFs + 1):2*N_DOFs);    % velocity        [m/s]
a0  = x0((2*N_DOFs + 1):3*N_DOFs);    % acceleration    [m/s2]
xi0 = x0((3*N_DOFs + 1):     end);    % BW displacement [m]

% assess the initial restoring forces
r0 = alpha_BW.*k.*diff([0; u0]) + ...    % spring force [kN] -- {1} Eq.54
     (1 - alpha_BW).*k.*xi0;
r0 = r0 - [r0(2:end);  0];               % DOFs force   [kN]

%% TIME-INTEGRATION:
% assign the initial conditions
u  = NaN(N_DOFs, Nt);      u(:, 1)  = u0;     % [m]
v  = NaN(N_DOFs, Nt);      v(:, 1)  = v0;     % [m/s]
a  = NaN(N_DOFs, Nt);      a(:, 1)  = a0;     % [m/s2]
xi = NaN(N_DOFs, Nt);      xi(:, 1) = xi0;    % [m]
r  = NaN(N_DOFs, Nt);      r(:, 1)  = r0;     % [kN]

% simulate the deterministic system to get data
for i = 2:Nt
    % for each time step,
    % 1) sample the force with a linear approximation
    f_tau = (1 - alpha_f)*f(:, i) + alpha_f*f(:, i - 1);    % [kN]

    % 2) initalize the displacements and velocities
    u(:, i) = u(:, i - 1);    % [m]
    v(:, i) = v(:, i - 1);    % [m/s]

    % 3) compute the equivalent stiffness matrix and force vector
    K_hat = M*((1 - alpha_m)/(beta*dt^2))     + ...    % [kN/m]
            C*((1 - alpha_f)*gamma/(beta*dt));

    f_hat = f_tau                                                      - ... % [kN]
            (  M*(1 - 1/(2*beta) + alpha_m/(2*beta))                   + ...
               C*(dt*(1 - alpha_f)*(1 - gamma/(2*beta))) )*a(:, i - 1) - ...
            (  C*(1 - gamma/beta + alpha_f*gamma/beta)                 - ...
               M*((1 - alpha_m)/(beta*dt))               )*v(:, i - 1) - ...
            (- M*((1 - alpha_m)/(beta*dt^2))                           - ...
               C*((1 - alpha_f)*gamma/(beta*dt))         )*u(:, i - 1) - ...
            alpha_f*r(:, i - 1);

    % 4) find the new nodal displacements with Newton-Raphson iteration
    for iter = 1:max_iter
        % 4.1) calculate the BW response
        [xi(:, i), r(:, i), K] = BW_model(xi(:, i - 1),         ...    % [m], 
                                          v(:, i - 1), v(:, i), ...    % [kN],
                                          u(:, i), k, alpha_BW, ...    % [kN/m]
                                          beta_BW, gamma_BW, n_BW, dt);

        % 4.2) assess the tangent stiffness and compute the residual force
        K_t   = K_hat + (1 - alpha_f)*K;                            % [kN/m]
        R_for = f_hat - (K_hat*u(:, i) + (1 - alpha_f)*r(:, i));    % [kN]

        % 4.3) if tolerance is not met, update the current displacements
        %      and velocities
        if norm(R_for) > tolerance
            duN     = K_t\R_for;         % increment    [m]
            u(:, i) = u(:, i) + duN;     % displacement [m]

            v(:, i) = (gamma/(beta*dt))*(u(:, i) - u(:, i - 1)) + ... % velocity
                      (1 - gamma/beta)*v(:, i - 1)              + ... % [m/s]
                      dt*(1 - gamma/(2*beta))*a(:, i - 1);
        else
            break
        end
    end

    % 5) compute the current accelerations
    a(:, i) = (1/(beta*dt^2))*(u(:, i) - u(:, i - 1)) - ...    % [m/s2]
              (1/(beta*dt))*v(:, i - 1)               - ...
              (1/(2*beta) - 1)*a(:, i - 1);
end

%% DATA CREATION:
% define the state and get the measurements
x = [u; v; a; xi];    % state        [m], [m/s], [m/s2], [m] -- {1} Page 3
y = H*x;              % measurements [m], [m/s], [m/s2]      -- {1} Eq.11

% add noise to the measurements
noise_mean = zeros(Ny, 1);                              % [m],  [m/s],   [m/s2]
noise_var  = diag(noise_ratio*var(y, 0, 2));            % [m2], [m2/s2], [m2/s4]
y          = y + mvnrnd(noise_mean, noise_var, Nt)';    % [m],  [m/s],   [m/s2]
end