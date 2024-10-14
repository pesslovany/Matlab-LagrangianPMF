function [estimatedStates,var] = particleFilterEstimationMult(modelChoose, souradniceGNSS, hB, vysky)


    % Initialize model
    model = initModel(modelChoose, souradniceGNSS, hB, vysky);

    % Define the number of particles
    numParticles = model.noPart;

    % Define state transition function
    stateTransitionFcn = @(x, u, model) model.F * x + u + sqrt(model.Q) * randn(size(x));

    % Define measurement likelihood function
    measurementLikelihoodFcn = @(x, z, model) ...
        mvnpdf(z', model.hfunct(x, 0, 0)', model.R);

    % Initialize the particle filter object
    pf = particleFilter(stateTransitionFcn, measurementLikelihoodFcn);

    pf.ResamplingMethod = 'systematic';
    pf.ResamplingPolicy.TriggerMethod = 'ratio';
    pf.ResamplingPolicy.SamplingInterval = 1;
    pf.ResamplingPolicy.MinEffectiveParticleRatio = 2/3;

    % Set number of particles
    initialize(pf, numParticles, model.meanX0, model.varX0);

    % Run the particle filter for the number of time steps
    estimatedStates = zeros(model.nx, length(model.z));  % store estimates

    for k = 1:length(model.z)
        % Predict the next state
        predict(pf, model.u(:, k), model);  

        % Correct using the measurement
        correct(pf, model.z(:, k), model);

        % Estimate the state (mean of the particle states)
        estimatedStates(:, k) = pf.getStateEstimate; % Now pf.State gives you the estimated state
        var(:,k) = diag(pf.StateCovariance);
    end

end

