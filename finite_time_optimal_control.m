function [x, u] = finite_time_optimal_control(A, B, X, U, X_f, Z, Q, R, P, N, x0, options, quad)
% by Lukas Brunke, lukas.brunke@mail.utoronto.ca
% Solves a constrained finite time optimal control problem

% Get problem dimensions
state_dim = size(A, 2);
input_dim = size(B, 2);

% Define optimization variables
x = sdpvar(state_dim, N + 1);
u = sdpvar(input_dim, N);

% Initialize cost
cost = 0;

% Initial constraint
% Determine if Z is a Polyhedron, otherwise use single initial constraint
polyhedron_string = 'Polyhedron';
if length(class(Z)) == length(polyhedron_string)
    if all(class(Z) == 'Polyhedron')
        Z = x0 + Z;
        constraints = Z.A * x(:, 1) <= Z.b;
    else
        constraints = x(:, 1) == x0;
    end
else
    constraints = x(:, 1) == x0;
end
for i = 1 : N
    % Create stage cost
    if quad % quadratic cost
        cost = cost + x(:, i)' * Q * x(:, i) + ...
               u(:, i)' * R * u(:, i);
    else % linear cost
        cost = cost + Q * sum(abs(x(:, i))) + ...
               R * sum(abs(u(:, i)));
    end
    
    % Create constraints from system dyanmics
    constraints = [constraints; x(:, i + 1) == A * x(:, i) + B * u(:, i)];
    % Set input and state constraints
    constraints = [constraints; 
                   X.A * x(:, i) <= X.b;
                   U.A * u(:, i) <= U.b;];                   
end

% Add terminal cost
if quad % quadratic cost
    cost = cost + x(:, N + 1)' * P * x(:, N + 1);
else % linear cost
    cost = cost + P * sum(abs(x(:, N + 1)));
end

% Add terminal constraint
constraints = [constraints; X_f.A * x(:, N + 1) <= X_f.b];

% Solve optimization problem
problem = optimize(constraints, cost, options);
objective = double(cost);
if problem.problem ~= 0
    double(x)
    error('No solution found.')
end

