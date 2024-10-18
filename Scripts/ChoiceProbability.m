% clc
% clear all
% close all

rng('shuffle');

expData = loadExpData('./Data/expData_g0.3_c_uambi_0.1.mat');

trialMatrix = expData.trialMatrix;
trialResponses = expData.trialResponses; % Note: this is logical matrix
neuronPrefOrientations = expData.preferredOrientation;

ambiguousStimuliOrientation = 90; % In degrees
neuronPreferredDirection = (neuronPrefOrientations > ambiguousStimuliOrientation) * (-1) + (neuronPrefOrientations <= ambiguousStimuliOrientation) * (1); % What about neurons which have 90 degree orientation
ambiguousStimTrialsIDx = find(trialMatrix(:, 2) == ambiguousStimuliOrientation);
ambiguousStimTrials = (trialMatrix(:, 2) == ambiguousStimuliOrientation);

nNeurons = length(neuronPrefOrientations);
neuronsCPs = zeros(nNeurons, 2); % Stores choice probability and its significance

% Choice probability is computed for each neuron
for neuronIDx=1:nNeurons
    % Based on this neuron's preferred orientation get the preferred trials
    % and null trials for this neurons.
    preferredTrialsIDx = ambiguousStimTrials & ( trialMatrix(:, 4) == neuronPreferredDirection(neuronIDx) );
    nullTrialsIDx = ambiguousStimTrials & ( trialMatrix(:, 4) ~= neuronPreferredDirection(neuronIDx) );
    preferredTrialsFR = sum(squeeze(trialResponses(preferredTrialsIDx, neuronIDx, :)), 2); % TODO: Normalize by total time to get firing rate
    nullTrialsFR = sum(squeeze(trialResponses(nullTrialsIDx, neuronIDx, :)), 2); % TODO: Normalize by total time to get firing rate
    
    % Compute choice probability based on firign rates distributions of
    % preferred vs null trials.
    firingRates = [preferredTrialsFR; nullTrialsFR];
    labels = [ones(size(preferredTrialsFR)); zeros(size(nullTrialsFR))];
    [~,~,~,AUC] = perfcurve(labels, firingRates, 1); % Compute the ROC curve and the AUC (Area Under the Curve)
    
    % To evaluate the statistical significance of the Choice Probability (CP) 
    % calculated for each neuron, a permutation test is performed. The procedure 
    % involves the following steps:
    % 
    % 1. Shuffle the trial labels (e.g., decision or stimulus condition) to 
    %    break any real association between the labels and the corresponding 
    %    neural firing rates.
    % 
    % 2. Recalculate the CP based on the firing rates of the neurons and these 
    %    shuffled labels. This randomization process disrupts any meaningful 
    %    relationships between neural activity and the actual trial outcomes.
    % 
    % 3. Repeat the shuffling and CP calculation multiple times (e.g., 1000 
    %    permutations) to generate a distribution of CP values under the 
    %    assumption that no true relationship exists between firing rates and 
    %    trial labels.
    % 
    % 4. Compare the actual CP value (calculated with the true labels) to this 
    %    distribution of CP values from the shuffled data to determine its 
    %    significance. This comparison allows for an assessment of whether the 
    %    observed CP is significantly different from chance.
    % 
    % The resulting p-value is the proportion of shuffled CP values that are 
    % equal to or greater than the actual CP, providing a measure of the 
    % statistical significance of the neuron's contribution to the choice process.
    numPermutations = 1000;
    permutedAUCs = zeros(numPermutations, 1);
    
    for i = 1:numPermutations
        permutedLabels = labels(randperm(length(labels)));
        [~,~,~,permutedAUCs(i)] = perfcurve(permutedLabels, firingRates, 1);
    end
    
    % Neuron's CP is significant is P-value is either less then 0.05 or
    % greater than 95.
    thresholdCP_95 = prctile(permutedAUCs, 95);
    thresholdCP_05 = prctile(permutedAUCs, 5);
    
    neuronsCPs(neuronIDx, 1) = AUC;
    neuronsCPs(neuronIDx, 2) = (AUC > thresholdCP_95) | (AUC < thresholdCP_05) ;
end

sigCPs = neuronsCPs(neuronsCPs(:,2)==1, :);
sigCPPrefOrientation = neuronPrefOrientations(neuronsCPs(:,2)==1);


% ----------------------------------
% Plot results
% ----------------------------------

figure
sgtitle('With noise correlation');

subplot(1, 2, 1)
h1 = histogram(neuronsCPs(:, 1), 20);
axis square;
binEdges = h1.BinEdges;
h1.FaceColor = 'none';
h1.EdgeColor = 'k';

hold on
h2 = histogram(sigCPs(:, 1), 'BinEdges', binEdges);
h2.FaceColor = 'k'; 
h2.EdgeColor = 'k';

legend('All CPs', 'Significant CPs')
legend boxoff; 
xlabel("Choice probability")
ylabel("Count")
xlim([0, 1])
ylim([0, 80])
title("CP distribution (all neurons)")
hold off

x = neuronPrefOrientations;
y = neuronsCPs(:,1);

p = polyfit(x, y, 30);  % You can change the degree (e.g., 1 for linear, 3 for cubic)
% Generate fitted values using the polynomial
xFit = linspace(min(x), max(x), 100);  % Create x values for plotting the fit
yFit = polyval(p, xFit);  % Evaluate the polynomial at the fitted x values

subplot(1, 2, 2)
hold on;
box off; axis square;
% xlim([70 110]); %xlim([45 135]);
ylim([0.3 0.9]); %ylim([0 1]);
s1 = scatter(x, y);
plot(xFit, yFit, '-b', 'LineWidth', 1);
xlabel("Preferred orientation (deg)")
ylabel("CP")
xline(90, '--r', 'LineWidth', 1);
s2 = scatter(sigCPPrefOrientation, sigCPs(:,1), 'filled', 'MarkerFaceAlpha', 0.5);
legend(s2, {'Significant CPs'});
legend boxoff; 
hold off;

%%
sigCP_70_80 = sigCPs(sigCPPrefOrientation > 70 & sigCPPrefOrientation < 80, 1);
sigCP_15_35 = sigCPs(sigCPPrefOrientation > 15 & sigCPPrefOrientation < 35, 1);

[p1, h1, stats1] = ranksum(sigCP_15_35, sigCP_70_80);

CP_70_80 = neuronsCPs(neuronPrefOrientations > 70 & neuronPrefOrientations < 80, 1);
CP_80_100 = neuronsCPs(neuronPrefOrientations > 80 & neuronPrefOrientations < 100, 1);

[p2, h2, stats2] = ranksum(CP_70_80, CP_80_100);


%% Correlation b/w previous trials and current ambiguous trials
rng('shuffle');

ambiguousTrialsIDX = find(trialMatrix(:, 2) == ambiguousStimuliOrientation);
prevTrialsIDx = ambiguousTrialsIDX - 1;
decisionAmbiguousTrials = trialMatrix(ambiguousTrialsIDX, 4);
decisionPrevTrials = trialMatrix(prevTrialsIDx, 4);

figure
subplot(2, 2, 1)
x = decisionAmbiguousTrials + 0.5*rand(size(decisionAmbiguousTrials));
y = decisionPrevTrials + 0.5*rand(size(decisionPrevTrials));
scatter(x, y)
xlabel("Decision of ambiguous trials")
ylabel("Decision of preceeding trials")
xlim([-1.3, 1.3])
ylim([-1.3, 1.3])

corr(x, y) 
% Gain 0.3 = 0.4262 Gain = -0.0728, Gain 0.6 = 0.8462, Gain 0.1=0.0788, 
% Gain=0.2 = 0.4267, Gain 0.4 = 0.6511, Gain 0.5 = 0.7398
% Gain -0.6= -0.6239, -0.5 = -0.8069, -0.4 = -0.74, -0.3 = -0.23, -0.2 =
% -0.45, -0.1 = -0.1808

gains = [-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
Corr = [-0.6239, -0.8069, -0.74, -0.23, -0.45, -0.1808, -0.0728, 0.0788, 0.4267, 0.4262, 0.6511, 0.7398, 0.8462];

subplot(2, 2, 2)
scatter(gains, Corr)
xlabel("Top-down gain value")
ylabel("Corr (ambiguous trials "+ newline +"with preceeding trials)")
xlim([-0.8, 0.8])