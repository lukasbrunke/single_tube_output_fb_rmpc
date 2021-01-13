function [X_f, nu] = terminal_constraint(X, U, A, B, K)
% by Lukas Brunke, lukas.brunke@mail.utoronto.ca
% According to Kouvaritakis and Cannon (2016), section 2.4 Theorem 2.3

state_dim = size(A, 2);
input_dim = size(B, 2);

% Rewrite constraints in joint form: F * x + G * u <= 1
[F, G, n_c] = to_joint_form(X, U);

% Define Phi
Phi = A + B * K;
Phi_k = eye(length(Phi));

F_GK = F + G * K;

% Initialize Phi^k, 
Phi_k = Phi_k * Phi;
mrpi = Polyhedron(F_GK * Phi_k, ones(n_c, 1));
mrpi.minHRep();
nu = 0;

while true
    nu = nu + 1;
    Phi_k = Phi_k * Phi;
    mrpi_next = Polyhedron(F_GK * Phi_k, ones(n_c, 1));
    mrpi_next.minHRep();
    mrpi_tmp = mrpi.intersect(mrpi_next);
    mrpi_tmp.minHRep();
    if mrpi_tmp.contains(mrpi) && mrpi.contains(mrpi_tmp)
        disp('Found mRPI')
        break
    else
        mrpi = mrpi_tmp;
    end
    
    if nu > 20
        disp('Reached maximum number of iterations. Did not find mRPI')
        break
    end
end

% Final nu
nu

Phi_k = eye(length(Phi));
num_constraints = length(F_GK);
V_T = zeros(length(F_GK) * nu, state_dim); 
for i = 0 : nu
    if i > 0
        Phi_k = Phi_k * Phi;
    end
    V_T(1 + i * num_constraints : (1 + i) * num_constraints, :) = F_GK * Phi_k;
end

% Calculate terminal set as intersection of all sets 
% (F + G * K) * Phi^i <= 1, for all i = 0, ..., nu
X_f = Polyhedron(V_T, ones(num_constraints * (nu + 1), 1));
X_f.minHRep();
end

