function [rmse, astd, annes] = calcRes(x, measMean, measVar, k)
%GBFEVALUATE  Evaluate RMSE, average std, and ANEES for GbF (any dimension)
%
%   [rmsePMF_col, astd, annesOut] = ...
%       gbfEvaluate(x, measMean2, measVar2, k)
%
%   INPUTS:
%     x           : true states       [nx × k]
%     measMean2   : filter means      [nx × k]
%     measVar2    : filter variances  [nx × nx × k]
%     k           : current time index (evaluate from 1 to k-1)
%
%   OUTPUTS:
%     rmsePMF_col : [nx × 1] RMSE vector
%     astd        : [nx × 1] sqrt of mean of diagonal entries of measVar2
%     annesOut    : scalar, accumulated ANEES

    nx = size(x,1);
    rmse = sqrt(mean((x(:,1:k-1) - measMean(:,1:k-1)).^2, 2));

    astd = zeros(nx,1);
    for i = 1:nx
        astd(i) = sqrt(mean(measVar(i,i,1:k-1)));
    end

    annes_PMF = 0;
    for indAn = 1:k-1
        diff = x(:,indAn) - measMean(:,indAn);
        invDiag = 1 ./ diag(measVar(:,:,indAn));
        annes_PMF = annes_PMF + (diff .* invDiag)' * diff;
    end
    annes = annes_PMF;
    
end
