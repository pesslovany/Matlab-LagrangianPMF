function A = setEdgesToZeros(A)
    % Set edges of PMD to zeros to ensure stability of the spectral
    % differentiation
    
    dims = ndims(A); % Get the number of dimensions of matrix A
    sizeA = size(A); % Get the size of matrix A
    
    % Loop over each dimension to set the edges to zero
    for dim = 1:dims
        % Create index vectors for all dimensions
        index = repmat({':'}, 1, dims);
    
        % Set the first edge along this dimension
        index{dim} = 1; % Update index to point to the first element
        A(index{:}) = 0; % Set the first edge to zero
    
        % Set the last edge along this dimension
        index{dim} = sizeA(dim); % Update index to point to the last element
        A(index{:}) = 0; % Set the last edge to zero
    end
end