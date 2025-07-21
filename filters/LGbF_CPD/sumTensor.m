function s = sumTensor(X)
%SUMTENSOR Sum of all elements in a ktensor (CPD)
%   s = sumTensor(X) returns the sum over all elements of the CPD tensor X,
%   without converting it to a full tensor.
%
%   Equivalent to sum(X(:)) for full tensors.

    assert(isa(X, 'ktensor'), 'Input must be a ktensor');

    R = length(X.lambda);
    N = ndims(X);
    s = 0;

    for r = 1:R
        prod_term = X.lambda(r);
        for n = 1:N
            prod_term = prod_term * sum(X.u{n}(:,r));
        end
        s = s + prod_term;
    end
    
end
