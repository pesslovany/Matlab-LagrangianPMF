function [cpd] = gbfDiffusionCPD(...
    predDensityCPD, Q, predGridDeltaCPD, nx, k, Rdef)
% Discrete diffusion by convolution with noise kernel

  % Grid dx
  dx_vec = predGridDeltaCPD(:,k+1);

  % Create kernel and convolve
  Uhat = cell(nx,1);
  for d = 1:nx
    sigma = sqrt(Q(d,d));       
    dx    = dx_vec(d);          
    sigma_i = sigma / dx;
    radius  = ceil(6 * sigma_i);
    x       = (-radius:radius) * dx;

    % Noise kernel
    G = exp(-x.^2/(2*sigma^2));
    G = G / (sum(G) * dx);

    % impose Dirichlet BC
    U = predDensityCPD.U{d};

    % 1D convolution along the first dim, same size
    Uhat{d} = convn(U, G', 'same');
  end

  % Deflate
  cpd = ktensor(predDensityCPD.lambda, Uhat);

  % Normalize
  Zsum       = sumTensor(cpd);
  cpd.lambda = cpd.lambda / (Zsum * prod(dx_vec));

end