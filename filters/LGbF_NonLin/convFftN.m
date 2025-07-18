function [pdfIn, kernel] = convFftN(pdfIn, kernel, Npa, kernelFFT)
% Calculates convolution by FFT using fftn and ifftn
% INPUTS:
% pdfIn     - pdf for convolution
% kernel    - kernel for convolution
% Npa       - number of points per axis (vector of length nx)
% nx        - dimension
% kernelFFT - is kernel already in frequency space?
% OUTPUTS:
% pdfIn     - convolution result

% Compute padded size
padSize = 2 * Npa - 1;

% Compute FFT of inputs
pdfIn = fftn(pdfIn, padSize);
if ~kernelFFT
    kernel = fftn(kernel, padSize);
end

% Perform convolution in frequency domain
pdfIn = pdfIn .* kernel;

% Back to state space
pdfIn = ifftn(pdfIn);

% Extract central part (remove padding)
subs = arrayfun(@(n) ceil((n-1)/2)+(1:n), Npa, 'UniformOutput', false);
pdfIn = real(pdfIn(subs{:}));

end