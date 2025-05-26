function delta_array = deltaToArray(delta_matrix)
    [J, ~, M] = size(delta_matrix);
    delta_array = [];
    for m = 1:M
        % Extract the J*J matrix at the m-th slice
        current_matrix = delta_matrix(:, :, m);
        
        % Get the linear indices of the lower triangular part, excluding the diagonal
        lower_triangle_indices = find(tril(true(J), -1));
        
        % Extract the lower triangular elements using the indices
        lower_triangle_vector = current_matrix(lower_triangle_indices);
        
        % Concatenate the lower triangular elements into the result array
        delta_array = [delta_array; lower_triangle_vector];
    end
end