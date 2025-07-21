function Z = hadamard_ktensor(A, B)
% Fast Hadamard product of two ktensors
%   Z = hadamard_ktensor(A, B) returns a ktensor Z = A .* B
%   (rank R = R_A * R_B) with all factor matrices built in one go.

    Ra = numel(A.lambda);
    Rb = numel(B.lambda);
    N  = ndims(A);

    %--- precompute index maps and new lambdas once ---
    [rIdx, sIdx] = ndgrid(1:Ra, 1:Rb);
    rIdx = rIdx(:);
    sIdx = sIdx(:);
    lambdaZ = A.lambda(rIdx) .* B.lambda(sIdx);

    factorsZ = cell(1, N);
    for n = 1:N
        Ua = A.u{n};   % size [I_n x Ra]
        Ub = B.u{n};   % size [I_n x Rb]

        % pick columns in one go, then elementwise multiply
        % Ua(:, rIdx) is [I_n x (Ra*Rb)], Ub(:, sIdx) same size
        factorsZ{n} = Ua(:, rIdx) .* Ub(:, sIdx);
    end

    Z = ktensor(lambdaZ, factorsZ);
    
end
