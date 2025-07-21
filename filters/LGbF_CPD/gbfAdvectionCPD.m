function [predDensityCPD, gridCPD, predGridDeltaCPD] = gbfAdvectionCPD(measPdfCPD, measMeanCPD, measVarCPD, ...
    gridCPD, F, Q, Npa, sFactor, predGridDeltaCPD, Rdef, k)
%GBFADVECTIONCPD calculates advection by point imaginary point movement and
%immediate interpolation back to rectangular axis aligne grid

% First two predictive moments by KF
diagAdvTrue = diag( F*diag(measVarCPD)*F' ) + diag(Q);     % 4×1
muPred      = F * measMeanCPD;                             % 4×1
stdPred     = sqrt(diagAdvTrue);

% Create axis aligned grid for these moments
axisp  = cell(4,1);
for d = 1:4
    span     = sFactor*stdPred(d);
    axisp{d} = linspace(muPred(d)-span, muPred(d)+span, Npa(d));
    predGridDeltaCPD(d,k+1) = axisp{d}(2) - axisp{d}(1);  % save Δx
end


% Build grid for each decomposable block of F
[ GX , GVx ] = ndgrid(axisp{1}, axisp{3});   % (x , ẋ) block
[ GY , GVy ] = ndgrid(axisp{2}, axisp{4});   % (y , ẏ) block

% The advection is solved as x_k+1 = F*x_k, this create a rotated grid
% which we cannot have, therefore we need to interpolate from x k+1 to
% aligned grid, however that is complicated, equivalent operation is to
% interpolated from x_k to F^-1 x_k+1, therefore form F^-1 x_k+1 grid
Xold  = GX  - GVx;
Vxold = GVx;
Yold  = GY  - GVy;
Vyold = GVy;

% Prepare emtpy loading vectors and lambdas
Unew       = {zeros(Npa(1),0), zeros(Npa(2),0), zeros(Npa(3),0), zeros(Npa(4),0)};
lambda_new = [];

uH       = measPdfCPD.U;
lambda   = measPdfCPD.lambda;

% Loop over loading vectors
for r = 1:ncomponents(measPdfCPD)
    % First decomposable slice, position velocity x interpolate
    u1s = interp1(gridCPD{1},uH{1}(:,r),Xold(:) ,'linear',0);
    u3s = interp1(gridCPD{3},uH{3}(:,r),Vxold(1,:),'linear',0);
    Zx  = reshape(u1s, Npa(1), Npa(3)) .* repmat(u3s,[Npa(1) 1]);          % Nx × Nvx

    % Decompose result to CPD loading vectors
    [Ux,Sx,Vx] = svd(Zx,'econ');
    sx     = diag(Sx);
    cumEx  = cumsum(sx.^2)/sum(sx.^2);
    Kx     = find(cumEx>=0.99999,1,'first');
    Ux_tr  = Ux(:,1:Kx)*sqrt(Sx(1:Kx,1:Kx));
    Vx_tr  = Vx(:,1:Kx)*sqrt(Sx(1:Kx,1:Kx));

    % Second decomposable slice, position velocity y interpolate
    u2s = interp1(gridCPD{2},uH{2}(:,r),Yold(:) ,'linear',0);
    u4s = interp1(gridCPD{4},uH{4}(:,r),Vyold(1,:),'linear',0);
    Zy  = reshape(u2s, Npa(2), Npa(4)) .* repmat(u4s,[Npa(2), 1]);          % Ny × Nvy

    % Decompose result to CPD loading vectors
    [Uy,Sy,Vy] = svd(Zy,'econ');
    sy     = diag(Sy);
    cumEy  = cumsum(sy.^2)/sum(sy.^2);
    Ky     = find(cumEy>=0.99999,1,'first');
    Uy_tr  = Uy(:,1:Ky)*sqrt(Sy(1:Ky,1:Ky));
    Vy_tr  = Vy(:,1:Ky)*sqrt(Sy(1:Ky,1:Ky));

    % Add all ranks together to form resulting advection solution on the
    % new grid
    for kx = 1:width(Ux_tr)
        for ky = 1:width(Uy_tr)
            Unew{1}(:,end+1) = Ux_tr(:,kx);   % x-axis
            Unew{3}(:,end+1) = Vx_tr(:,kx);   % ẋ-axis
            Unew{2}(:,end+1) = Uy_tr(:,ky);   % y-axis
            Unew{4}(:,end+1) = Vy_tr(:,ky);   % ẏ-axis
            lambda_new(end+1) = lambda(r);    %#ok<AGROW> % keep original weight
        end
    end
    % Not a great results when deflating each step
    % tempCP = ktensor(lambda_new', Unew);
    % if ncomponents(tempCP) > Rdef
    %     tempCP = cp_als(tempCP, Rdef);
    %     lambda_new = tempCP.lambda';
    %     for d = 1:4
    %         Unew{d} = tempCP.U{d};
    %     end
    % end
end

% New grid new weights
gridCPD         = axisp;                       
predDensityCPD  = ktensor(lambda_new',Unew);

% Deflate
if ncomponents(predDensityCPD) > Rdef
    predDensityCPD = cp_nmu(predDensityCPD, Rdef);
end

% Normalize
Zsum            = sumTensor(predDensityCPD);
predDensityCPD.lambda = predDensityCPD.lambda ...
    / (Zsum * prod(predGridDeltaCPD(:,k+1)));


end