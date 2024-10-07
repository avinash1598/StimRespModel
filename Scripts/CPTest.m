clc
clear all
close all

rng('shuffle');

expData = load('expData.mat');

trialMatrix = expData.expData.trialMatrix;
trialResponses = expData.expData.trialResponses;
neuronPrefOrientations = expData.expData.preferredOrientation;

trialStimOrientations = trialMatrix(:, 2);
stimVector = unique(trialStimOrientations);
countPerStim = length(trialStimOrientations) / length(stimVector);

ambiguousStimuliOrientation = 90; % In degrees
% For now directly set the preferred direction (-1 CCW, 1 CW)
neuronPreferredDirection = (neuronPrefOrientations > 90) * (-1) + (neuronPrefOrientations <= 90) * (1); % What about neurons which have 90 degree orientation
ambiguousStimTrialsIDx = find(trialMatrix(:, 2) == ambiguousStimuliOrientation);
ambiguousStimTrials = (trialMatrix(:, 2) == ambiguousStimuliOrientation);


neuronIDx = 7;
preferredTrialsIDx = find(ambiguousStimTrials & ( trialMatrix(:, 3) == neuronPreferredDirection(neuronIDx) ));
nullTrialsIDx = find(ambiguousStimTrials & ( trialMatrix(:, 3) ~= neuronPreferredDirection(neuronIDx) ));
% TODO: Normalize by total time to get firing rate
preferredTrialsFR = sum(squeeze(trialResponses(preferredTrialsIDx, neuronIDx, :)), 2);
nullTrialsFR = sum(squeeze(trialResponses(nullTrialsIDx, neuronIDx, :)), 2);

% Compute choice probability
firingRates = [preferredTrialsFR; nullTrialsFR];
labels = [ones(size(preferredTrialsFR)); zeros(size(nullTrialsFR))];
[~,~,~,AUC] = perfcurve(labels, firingRates, 1); % Compute the ROC curve and the AUC (Area Under the Curve)

% Display the AUC (Choice Probability)
fprintf('Area Under the ROC Curve (AUC): %.4f\n', AUC);

% Perform permutation test
numPermutations = 1000;
permutedAUCs = zeros(numPermutations, 1);

for i = 1:numPermutations
    permutedLabels = labels(randperm(length(labels)));
    [~,~,~,permutedAUCs(i)] = perfcurve(permutedLabels, firingRates, 1);
end

thresholdCP_95 = prctile(permutedAUCs, 95);
fprintf('Threshold CP: %.4f\n', thresholdCP_95);


% Optional: Plot the ROC curve
[Xroc, Yroc] = perfcurve(labels, firingRates, 1);
figure;
plot(Xroc, Yroc, '-');
hold on
x = linspace(0, 1, 100);
plot(x, x, '--')
xlabel('False Positive Rate (NULL trials)');
ylabel('True Positive Rate (PREF trials)');
title(sprintf('ROC Curve (AUC = %.4f)', AUC));
hold off


% Sample data (replace with your actual arrays)
array1 = preferredTrialsFR;
array2 = nullTrialsFR;

minVal = min([min(array1), min(array2)]);
maxVal = max([max(array1), max(array2)]);
edges = linspace(minVal, maxVal, 30);  % 30 bins
counts1 = histcounts(array1, edges);
counts2 = histcounts(array2, edges);

figure
bar(edges(1:end-1), counts1, 'FaceAlpha', 0.5, 'FaceColor', 'b', 'EdgeColor', 'b'); hold on;
bar(edges(1:end-1), -counts2, 'FaceAlpha', 0.5, 'FaceColor', 'r', 'EdgeColor', 'r');
xlabel('Data values');
ylabel('IPS');
title('Two Distributions with Common Bins');
legend('PREF', 'NULL');
grid on;
hold off;




% Single neuron preferred orientation computation
% nNeurons = length(neuronPrefOrientations);
% neuronPreferredDirection = zeros(1, nNeurons);

% sigModulatedNeuron = zeros(1, nNeurons);

% for neuronIDx=1:nNeurons
%     % neuronIDx = 1;
%     tuningFn = zeros(1, length(stimVector));
%     tuningFnPerTrial = zeros(length(stimVector), countPerStim);
% 
%     for i=1:length(stimVector)
%         stimOrientation = stimVector(i);
%         givenOrientationsTrialIDxes = find(trialStimOrientations == stimOrientation);
% 
%         neuronResponsesAtGivenOrientation = squeeze(trialResponses(givenOrientationsTrialIDxes, neuronIDx, :));
% 
%         avgFRAtGivenOrientationPerTrial = sum(neuronResponsesAtGivenOrientation, 2);
%         tuningFnPerTrial(i, :) = squeeze(avgFRAtGivenOrientationPerTrial);
% 
%         avgFRAtGivenOrientation = sum(sum(neuronResponsesAtGivenOrientation, 2)) / length(givenOrientationsTrialIDxes);
%         tuningFn(i) = avgFRAtGivenOrientation;
%     end
% 
%     % Perform one-way ANOVA for this neuron
%     [p, tbl, stats] = anova1(tuningFnPerTrial, [], 'off');  % 'off' suppresses the figure
% 
%     % [h, p] = swtest(tuningFnPerTrial); % Shapiro-Wilk test for normality
% 
%     % Check if p-value is significant (e.g., p < 0.05)
%     if p < 0.05
%         sigModulatedNeuron(neuronIDx) = p;
%         % disp('The neuron shows significant modulation by disparity.');
%     else
%         sigModulatedNeuron(neuronIDx) = p;
%         % disp('No significant modulation by disparity.');
%     end
% end

% neuronIDx = 3;
% tuningFn = zeros(1, length(stimVector));
% tuningFnPerTrial = zeros(length(stimVector), countPerStim);
% 
% for i=1:length(stimVector)
%     stimOrientation = stimVector(i);
%     givenOrientationsTrialIDxes = find(trialStimOrientations == stimOrientation);
% 
%     neuronResponsesAtGivenOrientation = squeeze(trialResponses(givenOrientationsTrialIDxes, neuronIDx, :));
% 
%     avgFRAtGivenOrientationPerTrial = sum(neuronResponsesAtGivenOrientation, 2);
%     tuningFnPerTrial(i, :) = squeeze(avgFRAtGivenOrientationPerTrial);
% 
%     avgFRAtGivenOrientation = sum(sum(neuronResponsesAtGivenOrientation, 2)) / length(givenOrientationsTrialIDxes);
%     tuningFn(i) = avgFRAtGivenOrientation;
% end
% 
% figure
% hold on
% scatter(stimVector, tuningFn);
% % for i=1:countPerStim
% %     plot(stimVector, tuningFnPerTrial(:, i));
% % end
% scatter(stimVector, sum(tuningFnPerTrial, 2)/countPerStim);
% xlabel("Orientations (deg)")
% ylabel("IPS")
% hold off
% 
% % Given that I have preferred vs null direction of neuron what to do next?
% % Choose range of ambiguous stimuli - based on fitted sigmoidal curve
% % Calculate choice eprobability of ambiguous stimuli
% 
% 
% % For CP just consider 90 for now