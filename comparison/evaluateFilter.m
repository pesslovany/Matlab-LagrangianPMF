function [rmse, astd11, astd22, astd33, annesOut] = evaluateFilter(x, estimatedState, covariance, method, k)
    % Function to evaluate a filtering method
    % Inputs:
    % x - true state or measurement
    % estimatedState - estimated state matrix (varies for each filter method)
    % covariance - variance matrix (varies for each filter method)
    % method - type of filter method (e.g., 'Sys', 'Mult', 'PMF', etc.)
    % k - current index for evaluation
    %
    % Outputs:
    % rmse - Root Mean Square Error
    % astd11, astd22, astd33 - Average standard deviations
    % annesOut - Annes output value

    % Calculate RMSE
    rmse = sqrt(mean((x(:, 1:k-1) - estimatedState(:, 1:k-1)).^2, 2));

    % Calculate average standard deviations
    if strcmp(method, 'PMF')
        astd11 = sqrt(mean(covariance(1, 1, 1:k-1)));
        astd22 = sqrt(mean(covariance(2, 2, 1:k-1)));
        astd33 = sqrt(mean(covariance(3, 3, 1:k-1)));
    else
        astd11 = sqrt(mean(covariance(1, 1:k-1)));
        astd22 = sqrt(mean(covariance(2, 1:k-1)));
        astd33 = sqrt(mean(covariance(3, 1:k-1)));
    end

    % Calculate Annes output
    annes = 0;
    for indAn = 1:k-1
        if strcmp(method, 'PMF')
            annes = annes + ((x(:, indAn) - estimatedState(:, indAn)) .* (1 ./ diag(covariance(:, :, indAn))))' * (x(:, indAn) - estimatedState(:, indAn));
        else
            annes = annes + ((x(:, indAn) - estimatedState(:, indAn)) .* (1 ./ (covariance(:, indAn))))' * (x(:, indAn) - estimatedState(:, indAn));
        end
    end
    annesOut = annes;
end
