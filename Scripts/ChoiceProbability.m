% clc
% clear all
% close all

rng('shuffle');

expData = loadExpData('./Data/expData_modgain_0.9_inv.mat');

trialMatrix = expData.trialMatrix;
trialResponses = expData.trialResponses;
neuronPrefOrientations = expData.preferredOrientation;

ambiguousStimuliOrientation = 90; % In degrees
neuronPreferredDirection = (neuronPrefOrientations > 90) * (-1) + (neuronPrefOrientations <= 90) * (1); % What about neurons which have 90 degree orientation
ambiguousStimTrialsIDx = find(trialMatrix(:, 2) == ambiguousStimuliOrientation);
ambiguousStimTrials = (trialMatrix(:, 2) == ambiguousStimuliOrientation);

nNeurons = length(neuronPrefOrientations);
neuronsCPs = zeros(nNeurons, 2); % Stores choice probability and its significance

for neuronIDx=1:nNeurons
    preferredTrialsIDx = ambiguousStimTrials & ( trialMatrix(:, 4) == neuronPreferredDirection(neuronIDx) );
    nullTrialsIDx = ambiguousStimTrials & ( trialMatrix(:, 4) ~= neuronPreferredDirection(neuronIDx) );
    preferredTrialsFR = sum(squeeze(trialResponses(preferredTrialsIDx, neuronIDx, :)), 2); % TODO: Normalize by total time to get firing rate
    nullTrialsFR = sum(squeeze(trialResponses(nullTrialsIDx, neuronIDx, :)), 2); % TODO: Normalize by total time to get firing rate
    
    % Compute choice probability
    firingRates = [preferredTrialsFR; nullTrialsFR];
    labels = [ones(size(preferredTrialsFR)); zeros(size(nullTrialsFR))];
    [~,~,~,AUC] = perfcurve(labels, firingRates, 1); % Compute the ROC curve and the AUC (Area Under the Curve)
    
    % Perform permutation test
    numPermutations = 1000;
    permutedAUCs = zeros(numPermutations, 1);
    
    for i = 1:numPermutations
        permutedLabels = labels(randperm(length(labels)));
        [~,~,~,permutedAUCs(i)] = perfcurve(permutedLabels, firingRates, 1);
    end
    
    thresholdCP_95 = prctile(permutedAUCs, 95);
    thresholdCP_05 = prctile(permutedAUCs, 5);
    
    neuronsCPs(neuronIDx, 1) = AUC;
    neuronsCPs(neuronIDx, 2) = (AUC > thresholdCP_95) | (AUC < thresholdCP_05) ;
end

sigCPs = neuronsCPs(neuronsCPs(:,2)==1, :);
sigCPPrefOrientation = neuronPrefOrientations(neuronsCPs(:,2)==1);

%%
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
ylim([0, 150])
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
xlim([0 180]); %xlim([45 135]);
ylim([0 1]); %ylim([0 1]);
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