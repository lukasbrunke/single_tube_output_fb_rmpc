function [F, G, n_c] = to_joint_form(X, U)
% by Lukas Brunke, lukas.brunke@mail.utoronto.ca
% Rewrite constraints in joint form: F * x + G * u <= 1
F = X.A;
f = X.b;
G = U.A;
g = U.b;
num_state_constraints = length(f);
num_input_constraints = length(g);
for i = 1 : num_state_constraints
    F(i, :) = 1 / f(i) * F(i, :);
end
for i = 1 : num_input_constraints
    G(i, :) = 1 / g(i) * G(i, :);
end

state_dim = size(F, 2);
input_dim = size(G, 2);

n_c = num_state_constraints + num_input_constraints;
new_F = zeros(n_c, state_dim);
new_G = zeros(n_c, input_dim);
new_F(1 : num_state_constraints, :) = F;
new_G(num_state_constraints + 1 : end, :) = G;

% Assign new F and G matrices
F = new_F;
G = new_G;
end

