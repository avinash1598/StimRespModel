% TIP: Keep all possible 1D variables to row vector
% Note: Time step should be set carefully. It should be pretty small
% relative to the spike rate of single neuron so that the probability of
% firing does not shoot up.
% TODO: add a check to make sure firing rate is not so big compared to the
% time window.

% ----------------------------------
% Params
% ----------------------------------
nNeurons = 50;       % Count of neurons
stimDuration = 1;    % Stimulus duration set to 1 seconds
varGain = 0.5;       % Variance in gain for modulated poisson process
timeStep = 0.001;    % 0.001s (1ms) - Step size of time bins used for binning stimulus duration 
% ntrials = 100;     % Offloaded to generateStim function. No of stimulus to simulate the model with

stimParam.startInterval = -pi;
stimParam.endInterval = pi;
stimParam.numStim = 20;
stimParam.countPerStim = 10;
ntrials = stimParam.numStim*stimParam.countPerStim;

neuronsPrefOrientation = zeros(1, nNeurons);
stimVector = generateStimVector(stimParam); % Stimulus vector containing orientation for each trial
timeBins = 0:timeStep:stimDuration; 
stimRespProfile = [0, normcdf(timeBins(2:end), .25, .05) .* 1./timeBins(2:end)];
gainVector = repmat(gamrnd(1./varGain, varGain, [ntrials, nNeurons]), [1, 1, length(timeBins)]); % for modulated poisson process - gain for each trial, neuron, time bin

% Neurons tuning parameters
tuningParams.d = zeros(1, nNeurons) + 1;       % Direction selectivity - set it to zero (no need for neuron to be directional selective).
tuningParams.alpha = zeros(1, nNeurons) + 2;   % Aspect ratio - controls sharpness.
tuningParams.b = zeros(1, nNeurons) + 2;       % Control sharpness of the neuron. Set it to 2 for these simulations.
tuningParams.q = zeros(1, nNeurons) + 1;       % Set it to some constant. Controls the amplitude of peak FR.
tuningParams.w = zeros(1, nNeurons) + 1;       % Doesn't matter what is val is becz untuned filter amp is zero.
tuningParams.UNTUNED_FILTER_AMPL = 0;          % Untuned filter not needed.
tuningParams.eps1 = 0;                         % Controls dynamic range.
tuningParams.beta = 1;                         % Controls dynamic range.
% Ask how to set this?
% how to do the normalization?
% Do i have to consider so many other variables out there?

% Temporary code
neuronsPrefOrientation(:) = -pi + (2 * pi) * rand(1, nNeurons);   % Randomly choose neurons preferred orientation from -pi to pi
% neuronsPrefOrientation = [-2.3880, -1.6928];

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

% STEP3:
% Trial by trial processing - for added complexity at later point
% Drive outpout of modulated poisson process
trialIDx = 1;
trlStimResponse = squeeze(stimResponse(trialIDx, :, :));
trlGainVector = squeeze(gainVector(trialIDx, :, :));

params = struct();
params.timeStep = timeStep;
params.timeBins = timeBins;
params.nNeurons = nNeurons;
[spikes, modStimResponse] = generateModulatedPoissonSpikes(trlStimResponse, ...
    trlGainVector, params);

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


figure
subplot(2, 2, 1)
hold on
title("Single trial response" + newline + "all neurons")
for nIDx = 1:nNeurons
    plot(squeeze(timeBins), squeeze(modStimResponse(nIDx, :)));
end
axis square
xlabel("Time")
ylabel("Firing rate")
hold off

subplot(2, 2, 2)
hold on
title("Single trial spikes" + newline + "all neurons")
imagesc(squeeze(spikes(:, :))), colormap(flipud('gray'))
box off, axis square, axis off
xlabel("Time")
ylabel("Firing rate")
hold off

subplot(2,2,3)
hold on
plot(decodingError)
xlabel("time")
ylabel("Estimated theta" + newline + "(MLE)")
hold off

% Sanity check - Mean of modStimRep and FR from spikes should be same
numIntervals = 50;
intervalSize = floor(length(timeBins) / numIntervals);  % 20 columns per interval
meanStimResp = zeros(nNeurons, numIntervals);
meanSpkRate = zeros(nNeurons, numIntervals);

for i = 1:numIntervals
    startCol = (i-1) * intervalSize + 1;
    endCol = min(i * intervalSize, length(timeBins));  % Handle the last interval

    % mean spk rate
    intervalData = sum(spikes(:, startCol:endCol), 2); 
    rate = intervalData / (timeStep*intervalSize);
    meanSpkRate(:, i) = rate;

    % mean stim response
    intervalData = mean(modStimResponse(:, startCol:endCol), 2); 
    meanStimResp(:, i) = intervalData;
end

subplot(2, 2, 4)
hold on
title("Mean spk rate")
for nIDx = 1:nNeurons
    plot(squeeze(meanSpkRate(nIDx, :)));
end
axis square
xlabel("Time bins")
ylabel("Firing rate")
hold off

% subplot(2,2,4)
% hold on
% x = linspace(1, 100, 100);
% plot(x, x, 'LineStyle','--')
% scatter(meanStimResp(:), meanSpkRate(:))
% xlabel("Mean (Mod stim resp)")
% ylabel("Mean spk rate")
% axis square
% hold off