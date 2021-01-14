%% Robust Output Feedback MPC: Lorenzetti and Pavone (2020)
% by Lukas Brunke, lukas.brunke@mail.utoronto.ca

% Implementation of Lorenzetti and Pavone: A Simple and Efficient Tube-Based 
% Robust Output Feedback Model Predictive Control Scheme, (2020).

% Example taken from Mayne et al. (2006)

%% Settings
% Make results in experiment repeatable:
experiment_1 = false; % select either true or false

if experiment_1
    % Define intial condition
    x0 = [-8.7; -1.5]; 
    rng(0,'philox');
else
    % Define intial condition
    x0 = [-8.7; -3.5];
    rng(0,'twister');
end

% Function to sample a random vertex from a convex set
random_vertex = @(X) X.V(randperm(length(X.V), 1), :)';

% linear or quadratic optimization problem
quadratic = true;

% Solver settings
if quadratic
    % solver = 'gurobi';
    solver = 'quadprog';
else
    % solver = 'gurobi';
    solver = 'linprog';
end
options = sdpsettings('verbose',0,'solver',solver);
options.OptimalityTolerance = 1e-15;
options.StepTolerance = 1e-15;
if strcmp('gurobi', solver)
    options.gurobi.IterationLimit = 1000;
    options.gurobi.TimeLimit = 5;
end

% Stopping criteria for determining mRPI sets
eps_rpi = 1e-5;
k_max = 50;

% Toggle maximum disturbances
maximum_disturbance = true;

%% Define system
% x(k + 1) = A * x(k) + B * u(k) + w(k)
%     y(k) = C * x(k) + v
A = [1.0, 1.0; 0.0, 1.0];
B = [1.0; 1.0];
C = [1.0, 1.0];

state_dim = size(A, 2);
input_dim = size(B, 2);
output_dim = size(C, 1);

%% Define constraints
% Input constraints
u_min = - 3.0;
u_max = 3.0;
U = Polyhedron([u_min; u_max]);
% plot(U)

% State constraints
x_1_min = -50.0;
x_1_max = 3.0;
x_2_min = -50.0;
x_2_max = 3.0;
X = Polyhedron([x_1_max, x_2_max; x_1_max, x_2_min; x_1_min, x_2_min; x_1_min, x_2_max]);
% plot(X)

% State disturbance constraints
w_min = - 0.1;
w_max = 0.1;
W = Polyhedron([w_max, w_max; w_max, w_min; w_min, w_min; w_min, w_max]);
% plot(W)

% Output disturbance constraints
v_min = - 0.05;
v_max = 0.05;
V = Polyhedron([v_min; v_max]);
% plot(V)

%% Define objectives
Q = eye(state_dim);
R = 0.01;

%% Determine control and observer gains
% Determine LQR state feedback controller
% Determine feedback gain and terminal cost matrix P from LQR for the
% unconstrained problem for the nominal system (without disturbances)
[K_f, P, ~] = dlqr(A, B, Q, R);
K_f = - K_f;

K = - [1, 1];
L = - [1; 1];

%% Closed-Loop system matrices
% Define error system matrix
A_k = A + B * K;
% Define observer system matrix
A_l = A + L * C;

%% Single Tube from Kögel and Findeisen (2017)
F = [A_l, zeros(state_dim);
     - L * C, A_k]; 
G = [eye(state_dim), L;
     zeros(state_dim), -L];
D_V = [W.V, V.V(1) * ones(length(W.V), 1);
       W.V, V.V(2) * ones(length(W.V), 1)];
D = Polyhedron(D_V);
D.minHRep();
GD = G * D;
GD.minHRep();

% Calculate mRPI set
if isnilpotent(F)
    disp('Calculating mRPI for nilpotent system matrix.')
    Z = rpi_nilpotent(F, GD);
else
    % According to Rakovic et al. (2005), Algorithm 1
    disp('Calculating mRPI using approximation in Rakovic et al. (2005).')
    Z = rpi(F, GD, eps_rpi, k_max);
end

% Assert that RPI set has been computed
Z_next = F * Z + GD;
Z_next.minHRep();
assert(Z.contains(Z_next) && Z_next.contains(Z))

% Tighten state and input constraints by Z
X_bar = X - [eye(state_dim), eye(state_dim)] * Z;
X_bar.minHRep();
U_bar = U - [zeros(1, state_dim), K] * Z;
U_bar.minHRep();

%% Determine terminal set as maximal roubst positive invariant set 
% According to Kouvaritakis and Cannon (2016), section 2.4 Theorem 2.3
X_f = terminal_constraint(X_bar, U_bar, A, B, K_f);

%% Define tube-based MPC
% Assumption in Lorenzetti: x_hat = x_bar, z_0 in Z_infty 
Z_sliced = Z.slice([3, 4], zeros(2, 1)); % slice Z at [e1, e2, 0, 0]^T
sampled_point = random_vertex(Z_sliced);
z0 = [sampled_point; zeros(2, 1)];

% Check if initial value is contained in Z for Proposition 2. in Lorenzetti
% and Pavone (2020)
assert(Z.contains(z0))

x_hat = x0 + sampled_point;
s0 = x_hat;

% Define the number of steps
num_steps = 50;

% Define MPC horizon
N = 7;

% Initialize open- and closed-loop trajectories for plotting
X_traj = zeros(state_dim, num_steps, N + 1);
X_traj_actual = zeros(state_dim, num_steps);
X_traj_hat = zeros(state_dim, num_steps);
U_traj = zeros(input_dim, num_steps, N);

% Run in closed-loop
for i = 1 : num_steps
    X_traj_actual(:, i) = x0;
    X_traj_hat(:, i) = x_hat;
    
    % Solve finite time optimal control problem
    % Note: Setting Z = zeros(state_dim, 1) avoids the optimization over an 
    % initial set at the first time step
    [x, u] = finite_time_optimal_control(A, B, X_bar, U_bar, X_f, zeros(state_dim, 1), Q, R, P, N, s0, options, quadratic);
    
    % Determine optimal control input
    u_opt = K * (x_hat - double(x(:, 1))) + double(u(:, 1));
    
    % Sample disturbance from set W and set V
    if maximum_disturbance
        % selecte disturbance from random vertices of the sets W and V
        w = random_vertex(W);
        v = random_vertex(V);
    else
        w = rejection_sampling(W);
        v = rejection_sampling(V);
    end
    
    % Set initial condition for next iteration by simulating the system
    y = C * x0 + v;
    x0 = A * x0 + B * u_opt + w;
    s0 = double(x(:, 2));
    
    % Observer
    y_hat = C * x_hat;
    x_hat = A * x_hat + B * u_opt + L * (y_hat - y);
    
    % Save open-loop trajectories
    X_traj(:, i, :) = double(x);
    U_traj(:, i, :) = double(u);
end

%% Plot Results
% Plot state space trajectories
figure
% Plot state constraints
plot(X, 'color','white')
hold on

plot(X_f + [eye(state_dim), eye(state_dim)] * Z, 'color', 'lightgray')
plot(X_f, 'color', 'gray')

% Plot tube over the nominal trajectory, the estimated trajectory and true
% trajectory
for i = 1 : num_steps
    Z_i = [X_traj(1, i, 1); X_traj(2, i, 1)] + [eye(state_dim), eye(state_dim)] * Z;
    plot(Z_i, 'color', 'g')
    if ~Z_i.contains([X_traj_actual(1, i); X_traj_actual(2, i)])
        i
        disp('True state is not contained inside the tube.')
    end
end

% Plot nominal closed-loop trajectory
plot(X_traj(1, :, 1), X_traj(2, :, 1), 'ks-')
% Plot actual closed-loop trajectory of the system with disturbance
plot(X_traj_actual(1, :), X_traj_actual(2, :), 'bo--')
% Plot closed-loop trajectory of the estimated state
plot(X_traj_hat(1, :), X_traj_hat(2, :), 'rx:')

% Set figure limits
xlim([-22, 5])
ylim([-11, 5])

% Set labels and title
xlabel("$x_1$", 'Interpreter','latex')
ylabel("$x_2$", 'Interpreter','latex')
title("Robust Output Feedback MPC")

% Set legend
f=get(gca,'Children');
lgd = legend([f(1), f(2), f(3), f(3 + num_steps), ...
              f(4 + num_steps), f(5 + num_steps), f(6 + num_steps)], ...
             '$\hat{x}$','$x$','$\bar{x}$',...
             '$\bar{x} \oplus (\vec{I} ~ \vec{I}) \mathcal{Z}$',...
             '$\mathcal{X}_f$',...
             '$\mathcal{X}_f \oplus (\vec{I} ~ \vec{I}) \mathcal{Z}$',...
             '$\mathcal{X}$',...
             'Interpreter', 'latex', 'Location', 'southwest');
lgd.FontSize = 16;

% Plot trajectories over time
figure
hold on
plot(1: num_steps, 3 * ones(num_steps, 1), 'k-')
plot(1 : num_steps, X_traj_actual(1, :), 'r-')
plot(1 : num_steps, X_traj_actual(2, :), 'b-')
hold off

% Set labels and title
xlabel("$t$", 'Interpreter','latex')
ylabel("$x_1,~ x_2$", 'Interpreter','latex')
title("Robust Output Feedback MPC")

% Set legend
f=get(gca,'Children');
lgd = legend([f(1), f(2), f(3)], ...
             '$x_1$','$x_2$','$\mathcal{X}$',...
             'Interpreter', 'latex', 'Location', 'southeast');
lgd.FontSize = 16;