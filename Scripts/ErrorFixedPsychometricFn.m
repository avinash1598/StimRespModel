% clear all
% close all
% clc

rng('shuffle'); % Keep this to ensure that different rnadom numbers are generated everytime

% TIP: Keep all possible 1D variables to row vector
% Note: Time step should be set carefully. It should be pretty small
% relative to the spike rate of single neuron so that the probability of
% firing does not shoot up.
% TODO: add a check to make sure firing rate is not so big compared to the
% time window.

% ----------------------------------
% Params
% ----------------------------------
nNeurons = 100;      % Count of neurons
stimDuration = 1;    % Stimulus duration set to 1 seconds
varGain = 0.5;       % Variance in gain for modulated poisson process
timeStep = 0.001;    % 0.001s (1ms) - Step size of time bins used for binning stimulus duration 
% ntrials = 100;     % Offloaded to generateStim function. No of stimulus to simulate the model with

stimParam.startInterval = deg2rad(90-8); 
stimParam.endInterval = deg2rad(90+8); 
stimParam.numStim = 21;
stimParam.countPerStim = 200;
ntrials = stimParam.numStim*stimParam.countPerStim;

stimVector = repelem(linspace(stimParam.startInterval, stimParam.endInterval, ...
    stimParam.numStim), stimParam.countPerStim); % Stimulus vector containing orientation for each trial

neuronsPrefOrientation = zeros(1, nNeurons);
timeBins = 0:timeStep:stimDuration; 
stimRespProfile = [0, normcdf(timeBins(2:end), .25, .05) .* 1./timeBins(2:end)];
gainVector = repmat(gamrnd(1./varGain, varGain, [ntrials, nNeurons]), [1, 1, length(timeBins)]); % for modulated poisson process - gain for each trial, neuron, time bin

% Neurons tuning parameters
tuningParams.d = zeros(1, nNeurons) + 1;       % (fixed) Direction selectivity - set it to zero (no need for neuron to be directional selective).
tuningParams.alpha = zeros(1, nNeurons) + 2;   % (fixed) Aspect ratio - controls sharpness. Keep this fixed. Reducing the value makes the changes very rapid towards the end which we probably don't want.
tuningParams.b = zeros(1, nNeurons) + 2;       % (fixed maybe/variable - (0.5, some max - 3, 4 ...)) Control this - Control sharpness + range of the neuron. Set it to 2 for these simulations
tuningParams.q = zeros(1, nNeurons) + 1;       % (variable) Set it to some constant. Controls the sharpness and amplitude of peak FR.
tuningParams.w = zeros(1, nNeurons) + 1;       % (fixed) Doesn't matter what is val is becz untuned filter amp is zero.
tuningParams.UNTUNED_FILTER_AMPL = 0;          % (fixed) Untuned filter not needed.
tuningParams.eps1 = zeros(1, nNeurons);        % (variable) Controls dynamic range.
tuningParams.beta = zeros(1, nNeurons) + 1;    % (variable) Controls dynamic range.

% Overriding variable tuning parameters
% tuningParams.b(:) = 0.5 + 10*rand(1, nNeurons);        % Maybe change non linearity instead. This is derivative order so ideally should be integer. Range: 0.5-4
tuningParams.q(:) = 8*rand(1, nNeurons);                 % Range: 0 - 20 ips, gamrnd(2, 2, [1, nNeurons]);
tuningParams.beta(:) = lognrnd(2.5, 0.5, 1, nNeurons);     % lognrnd(mu, sigma, 1, nneurons)
tuningParams.eps1(:) = lognrnd(1, 0.8, 1, nNeurons);     % Range: 0 - 10 ips, Might not be correct initilization

% Temporary code
neuronsPrefOrientation(:) = pi * rand(1, nNeurons);                 % Randomly choose neurons preferred orientation from -pi to pi
% neuronsPrefOrientation(:) = -pi + (2 * pi) * rand(1, nNeurons);   % Randomly choose neurons preferred orientation from -pi to pi

% ----------------------------------
% Computing stimulus response begins
% ----------------------------------

% STEP1: 
firingRates = orientationTunedFiringRate(stimVector, ...
    neuronsPrefOrientation, tuningParams);

% STEP2: 
% Find the time response of all these neurons to input stimuli
% This firing rate matrix will be multipled with some time response
% function to get the final stimulus response of each neuron.
% Perform elementwise multipliction between stimulus dependent firing rates
% and stimulus response profile.
stimResponse = firingRates.*reshape(stimRespProfile, 1, 1, length(stimRespProfile));

resultMat = zeros(stimParam.numStim, stimParam.countPerStim);
neuronSpikeResponses = zeros(ntrials, nNeurons, length(timeBins));
estimatedMLEs = zeros(stimParam.numStim, stimParam.countPerStim);

for i = 1:stimParam.numStim
    disp(i)
    for j=1:stimParam.countPerStim
        trialIDx = (i-1)*stimParam.countPerStim + j;
        
        % STEP3:
        % Trial by trial processing - for added complexity at later point
        % Drive outpout of modulated poisson process
        trlStimResponse = squeeze(stimResponse(trialIDx, :, :));
        trlGainVector = squeeze(gainVector(trialIDx, :, :));
        
        params = struct();
        params.timeStep = timeStep;
        params.timeBins = timeBins;
        params.nNeurons = nNeurons;
        [spikes, modStimResponse] = generateModulatedPoissonSpikes(trlStimResponse, ...
            trlGainVector, params);

        neuronSpikeResponses(trialIDx, :, :) = squeeze(spikes); % Store spikes in matrix to save later

        % STEP4:
        % Decode orientation based on modulated spikes
        params = struct();
        params.timeStep = timeStep;
        params.timeBins = timeBins;
        params.nNeurons = nNeurons;
        params.stimDuration = stimDuration;
        
        thetaMLE = decodeOrientationFromSpikes(spikes, ...
            neuronsPrefOrientation, params, tuningParams);
        decodingError = thetaMLE - stimVector(trialIDx);
        
        estimatedMLEs(i, j) = thetaMLE(end);

        decision = (thetaMLE(end) > pi/2)*(-1) + (thetaMLE(end) <= pi/2)*(1); % Maybe this should be generic i.e. don't hardode pi/2
        % actualVal = (stimVector(trialIDx) > pi/2)*(1) + (stimVector(trialIDx) <= pi/2)*(-1);

        resultMat(i, j) = decision; %(decision == 1);
    end
end

% Orientation selectivity index for each neuron
thetas_ = linspace(0, pi, 180);
tuningFnFR_ = orientationTunedFiringRate(thetas_, neuronsPrefOrientation, tuningParams)';
OSI = abs(tuningFnFR_*(exp(1i*(2*thetas_))'))./sum(abs(tuningFnFR_), 2);

% ----------------------------------
% Store results/responses
% ----------------------------------
trialMatrix = zeros(ntrials, 3); % trial no, stim orientation, estimated decision

for i = 1:stimParam.numStim
    for j = 1:stimParam.countPerStim
        trial_idx = j + (i - 1)*stimParam.countPerStim;

        trialMatrix(trial_idx, 1) = trial_idx;                          % Trial index
        trialMatrix(trial_idx, 2) = rad2deg(stimVector(trial_idx));     % Stimulus orientation
        trialMatrix(trial_idx, 3) = resultMat(i, j);                    % Decision 
    end
end

expData.trialMatrix = trialMatrix;
expData.columnDescriptions.trialMatrix = {'Column 1: Trial Number', 'Column 2: Stimulus orientation (in degree)', 'Column 3: Decision (1: CW, -1: CCW)'};
expData.preferredOrientation = rad2deg(neuronsPrefOrientation);
expData.trialResponses = neuronSpikeResponses;
expData.columnDescriptions.trialResponses = {'Spike responses for ntrials x nNeurons x nTimeBins'};
expData.timeBins = timeBins; % TODO: store time bin information

% save('expData.mat', 'expData', '-v7.3');

% ----------------------------------
% Plot results
% ----------------------------------

% figure
% subplot(2, 2, 1)
% hold on
% histogram(OSI, 10, 'Normalization', 'probability')
% xlabel('Orientation Selectivity')
% ylabel('Count')
% hold off
% 
% % Plot all tuning curves
% subplot(2, 2, 2)
% hold on
% for i = 1:nNeurons
%     plot(rad2deg(thetas_), tuningFnFR_)
% end
% xlabel('Orientation (rad)')
% ylabel('ips')
% hold off


figure
subplot(2, 2, 1)
hold on
title("Single trial response" + newline + "all neurons")
for nIDx = 1:nNeurons
    plot(squeeze(timeBins), squeeze(modStimResponse(nIDx, :)));
end
axis square
xlabel("Time (s)")
ylabel("IPS")
hold off

subplot(2, 2, 2)
hold on
title("Single trial spikes" + newline + "all neurons")
imagesc(squeeze(spikes(:, :))), colormap(flipud('gray'))
axis square
box off, axis off
xlabel("Time (s)")
ylabel("Spikes")
hold off

subplot(2,2,3)
hold on
title("Single trial")
timePts = linspace(0, stimDuration, length(decodingError));
plot(timePts, decodingError)
xlabel("Time (s)")
ylabel("Decoding error")
hold off

uniqStim = unique(stimVector);
psychometricData = zeros(1, length(uniqStim));

for i=1:stimParam.numStim
    stimOrientation = uniqStim(i);
    givenOrientationTrialIDxes = find(stimVector == stimOrientation);
    decision = resultMat(i, :);
    
    percent_CW = length(find(decision == -1)) / length(decision);
    psychometricData(i) = percent_CW;
end

x = rad2deg(uniqStim);
y = psychometricData;

subplot(2,2,4)
scatter(x, y, 'DisplayName', 'Data')
% scatter(y)
xlabel("Orientation (deg)")
ylabel("% CCW")
% legend;
title('Psychometric Function Fit');
% ylim([0, 1.2])
hold off

% stimVals = unique(stimVector);
% proportion = sum(resultMat' == 1) ./ size(resultMat', 1);
% subplot(2,2,4)
% % hold on
% x = rad2deg(stimVals);
% y = proportion;
% % plot(x, y, '-o')
% % xlabel("Orientation (deg)")
% % ylabel("Percentage > pi/2")
% % hold off
% 
% sigmoid = @(a, b, x) 1 ./ (1 + exp(-(a * (x - b))));
% initialGuess = [1, mean(x)]; % Initial guess for a and b
% params = nlinfit(x, y, sigmoid, initialGuess);
% a = params(1); % Slope
% b = params(2); % Inflection point (midpoint)
% xFit = linspace(min(x), max(x), 100); % Create fine x values for plotting the fitted curve
% yFit = sigmoid(a, b, xFit);           % Compute the y values based on the fitted parameters
% 
% % Plot the original data and the fitted sigmoid
% scatter(x, y, 'o', 'DisplayName', 'Data');  % Original data
% hold on;
% plot(xFit, yFit, '-', 'DisplayName', 'Fitted Sigmoid'); % Fitted sigmoid curve
% xlabel('Stimulus (Degrees)');
% ylabel('Proportion');
% legend;
% title('Sigmoid Fit to Data');
% 
% % Sanity check - Mean of modStimRep and FR from spikes should be same
% numIntervals = 50;
% intervalSize = floor(length(timeBins) / numIntervals);  % 20 columns per interval
% meanStimResp = zeros(nNeurons, numIntervals);
% meanSpkRate = zeros(nNeurons, numIntervals);
% 
% for i = 1:numIntervals
%     startCol = (i-1) * intervalSize + 1;
%     endCol = min(i * intervalSize, length(timeBins));  % Handle the last interval
% 
%     % mean spk rate
%     intervalData = sum(spikes(:, startCol:endCol), 2); 
%     rate = intervalData / (timeStep*intervalSize);
%     meanSpkRate(:, i) = rate;
% 
%     % mean stim response
%     intervalData = mean(modStimResponse(:, startCol:endCol), 2); 
%     meanStimResp(:, i) = intervalData;
% end
% 
% % subplot(2, 2, 4)
% % hold on
% % title("Mean spk rate")
% % for nIDx = 1:nNeurons
% %     plot(squeeze(meanSpkRate(nIDx, :)));
% % end
% % axis square
% % xlabel("Time bins")
% % ylabel("Firing rate")
% % hold off
% 
% % subplot(2,2,4)
% % hold on
% % x = linspace(1, 100, 100);
% % plot(x, x, 'LineStyle','--')
% % scatter(meanStimResp(:), meanSpkRate(:))
% % xlabel("Mean (Mod stim resp)")
% % ylabel("Mean spk rate")
% % axis square
% % hold off
% 




function [thetaMLE] = decodeOrientationFromSpikes(spikes, neuronsPrefOrientation, params, tuningParams)
    timeBins = params.timeBins;
    nNeurons = params.nNeurons;    
    stimDuration = params.stimDuration;

    % Decoding params
    tauStep = 0.02; % Increment decoding window by 20 ms
    numIntervals = floor(stimDuration / tauStep);
    intervalSize = floor(length(timeBins) / numIntervals); 

    aggregateSpikeCounts = zeros(nNeurons, numIntervals);
    taus = zeros(1, numIntervals);
    thetaMLE = zeros(1, numIntervals);
    % avg_pdf = zeros(1, size(decThetas, 2));
    
    % Bayesian decoding on the firing maps to decode stimulus direction
    % Decode for each trial separately first
    decThetas = linspace(-pi, pi, 360);
    tuningFnFR = orientationTunedFiringRate(decThetas, neuronsPrefOrientation, tuningParams);
 
    % Do I need to consider the time dependent response as well! Let's see if
    % decoding works just with tuning function without any knowledge of time
    % dependent response (response profile is same for all so should not make much difference)
    % Bin spike counts into coarser intervals
   
    % Get aggregate spike count of each neuron in the given decoding window
    for i = 1:numIntervals
        % Increasing temporal windows
        startCol = 1; %(i-1) * intervalSize + 1; %1
        endCol = min(i * intervalSize, length(timeBins));  % Handle the last interval
        taus(i) = i*tauStep; % TODO: set tau as well
    
        intervalData = sum(spikes(:, startCol:endCol), 2);
        aggregateSpikeCounts(:, i) = intervalData;
    end
    
    % iterate for spike count in each time bin
    for t=1:size(aggregateSpikeCounts, 2)
        n = aggregateSpikeCounts(:, t);
        f_theta = tuningFnFR';
        
        % Probability desnity function estimation for theta distribution
        temp = taus(t)*f_theta + eps; % Adding a very small noise to avoid log going to NaN
        log_pdf = sum(n.*log(temp), 1) - sum(temp, 1) - sum(log(factorial(n)), 1); % Fixed error. Bingo! now decoding works! Factorial is pushing the value really high
        % pdf = exp(log_pdf); % Just don't convert to exp to avoid values from blowing up! But averaging pdf won't work with log pdfs
        % normalized_pdf = pdf / sum(pdf);
        
        [max_val, idx_max] = max(log_pdf); %max(normalized_pdf);
        thetaMLE(t) =  decThetas(idx_max);
    end
end

function [spikes, modStimResponse] = generateModulatedPoissonSpikes(trlStimFRResponse, trlGainVector, params)
    timeStep = params.timeStep;
    timeBins = params.timeBins;
    nNeurons = params.nNeurons;
    
    modStimResponse = trlGainVector.*trlStimFRResponse; 
    % Note: this way of calculating probability of spiking in each time bin
    % relies on the assumption that there exist a very small time bin dt within
    % which no two spikes can be fired. Hence, firing rate of the neurons
    % should not be very great compared to the size of the time bin.
    probSpk = modStimResponse.*timeStep;
    rndNumbers = rand(nNeurons, length(timeBins)); % Generate for each trial and each time bin
    spikes = rndNumbers < probSpk;
end

function [FR] = getUnormalizedTunedResponse(orientations, neuronsPrefOrientations, params)
    % Extract parameters from the structure
    d = params.d;
    alpha = params.alpha;
    b = params.b;
    q = params.q;
    w = params.w;
    UNTUNED_FILTER_AMPL = params.UNTUNED_FILTER_AMPL;
    
    % TODO: Add a check to make sure all these parameters are provided

    % Mean firing rates profiles of all the neurons from 
    % the given input stimuli.
    T1_ = orientations' - neuronsPrefOrientations;
    T2_ = 1 + 0.5*( sign(cos(T1_)) - 1).*d;
    T3_ = cos(T1_).*exp(-0.5*(cos(T1_).^2).*(1 - alpha.^2) ).^b;

    R_TUNED = T2_.*T3_; % Tuned linear filter response
    R_UNTUNED = w*UNTUNED_FILTER_AMPL; % Untuned filter response
    COMBINED_CHANNEL_OUTPUT = R_TUNED - R_UNTUNED;
    % TODO: Normalization leave it for later: Ask Robbe
    % Apply rectifier to the linear filter output
    COMBINED_CHANNEL_OUTPUT(COMBINED_CHANNEL_OUTPUT < 0) = 0;
    
    FR = COMBINED_CHANNEL_OUTPUT;
end

function [FR] = orientationTunedFiringRate(orientations, neuronsPrefOrientations, params)
    % TUNING_FN_FR: Tuning function of all the neurons used for
    % normalization
    % Extract parameters from the structure
    q = params.q;
    eps1 = params.eps1;
    beta = params.beta;  

    unnormalizedResp = getUnormalizedTunedResponse(orientations, neuronsPrefOrientations, params);
    
    % Normalize combined chanel output between zero and one for each neuron
    % Once the tuning curves are normalized, adjust the amplitudes to set
    % baseline and peak firing rates of individual neurons.
    % COMBINED_CHANNEL_OUTPUT : ntrials * nNeurons
    % Note: For normalization, unnormalized tuning function is required
    tuningFnOrientations = linspace(0, pi, 180);
    unormalizedTuningFn = getUnormalizedTunedResponse(tuningFnOrientations, neuronsPrefOrientations, params);
    minVals = min(unormalizedTuningFn);  % Minimum value of each column
    maxVals = max(unormalizedTuningFn);  % Maximum value of each column
    combinedChannelOp = (unnormalizedResp - minVals) ./ (maxVals - minVals);

    % Compute firign rate using power law non-linearity
    FR = eps1 + beta.*(combinedChannelOp.^q);
end