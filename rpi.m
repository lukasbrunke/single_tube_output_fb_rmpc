function [F_alpha_s, alpha, s] = rpi(A, W, eps, s_max)
% Compute the robust positively invariant set from Rakovic et al. (2005).

% Set the precision
% eps = 1e-3;

s = 0;
% s_max = 10;

state_dim = size(A, 2);

% Create standard basis vectors
e = eye(state_dim);

% Collect all the powers of A
A_pow_all = {};
A_pow_s = A^s;
A_pow_all{1} = A_pow_s;

% Get the disturbance constraints
f = W.A';
g = W.b;

% Initialize sum of supports
positive_support = zeros(state_dim, 1);
negative_support = zeros(state_dim, 1);

while 1
    s = s + 1;
    % Calculate next power of A
    A_pow_s = A_pow_s * A;
    A_pow_all{s + 1} = A_pow_s;
    
    % Determine alpha_o(s) from eq. (11) from Rakovic et al. (2005).
    alpha_o_all = zeros(length(g), 1);
    for i = 1 : length(g)
        alpha_o_all(i) = support_function(A_pow_s' * f(:, i), W) / g(i);
    end
    alpha = max(alpha_o_all);
    
    % Determine sum of supports from eq. (13) from Rakovic et al. (2005).
    for j = 1 : state_dim
       positive_support(j) = positive_support(j) + ...
                             support_function(A_pow_all{s}(j, :)', W);
       negative_support(j) = negative_support(j) + ...
                             support_function(- A_pow_all{s}(j, :)', W);
    end
    
    % Determine M(s) from eq. (13) from Rakovic et al. (2005).
    M = max([positive_support; negative_support]);

    % Precision criterion from eq. (14) from Rakovic et al. (2005)
    if s >= 2  % make sure it runs for at least two iterations 
        if alpha <= eps / (eps + M) 
            disp('RPI set computed')
            break
        end
    end
    
    % Maximum number of iterations criterion
    if s >= s_max
        disp('Maximum number of iterations reached')
        break
    end
end

% Print alpha and s
alpha
s

% Compute F_s and scale it to give F(alpha, s)
W.minHRep();
F_s = A_pow_all{1} * W;
F_s.minHRep();
% Compute the Minkowski sum for F_s from eq. (2) from Rakovic et al. (2005).
for i = 2 : s
    F_s = F_s + A_pow_all{i} * W;
    F_s.minHRep();
end

% Scale the set F_s with the factor from eq. (5) from Rakovic et al. (2005)
F_alpha_s = 1 / (1 - alpha) * F_s;

%Plot F(alpha, s)
% plot(F_alpha_s)
% xlabel("x_1")
% ylabel("x_2")
% title(["F(\alpha, s) \subseteq F_{\infty}",
%        join(["\alpha=", num2str(alpha), "; s=", num2str(s)])])
end

