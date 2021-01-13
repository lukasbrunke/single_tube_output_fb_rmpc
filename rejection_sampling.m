function [selected_point] = rejection_sampling(S)
    % by Lukas Brunke, lukas.brunke@mail.utoronto.ca
    % Sample a random point from a polyhedron using rejection sampling

    outer = S.outerApprox();
    num_vertices = size(outer.V, 1);
    
    selected_point = zeros(1, size(outer.V, 2));
    first_iteration = true;
    
    % Sample random points until point is contained in the set S
    while S.contains(selected_point') == 0 || first_iteration
        rand_values = rand(num_vertices, 1);
        percentages = rand_values / sum(rand_values);
        selected_point = percentages' * outer.V;
        if first_iteration
            first_iteration = false;
        end
    end
    selected_point = selected_point'; 
end

