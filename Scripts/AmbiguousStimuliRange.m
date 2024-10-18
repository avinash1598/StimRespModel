% clc
% clear all
% close all

expData = loadExpData('./Data/expData_contrast_0.3.mat');

trialMatrix = expData.trialMatrix;
trialResponses = expData.trialResponses;
neuronPrefOrientations = expData.preferredOrientation;

stimVector = trialMatrix(:, 2);
noisyStimVector = trialMatrix(:, 3);

unqStimOrientations = unique(stimVector)';
countPerStim = length(noisyStimVector) / length(unqStimOrientations);

psychometricData = zeros(1, length(unqStimOrientations));

for i=1:length(unqStimOrientations)
    stimOrientation = unqStimOrientations(i);
    givenOrientationTrialIDxes = find(stimVector == stimOrientation);
    decision = trialMatrix(givenOrientationTrialIDxes, 4);
    
    percent_CCW = length(find(decision == -1)) / length(decision);
    psychometricData(i) = percent_CCW;
end

x = unqStimOrientations';
y = psychometricData';

figure
scatter(x, y, 'DisplayName', 'Data')
xlabel("Orientation (deg)")
ylabel("% CCW")
title('Psychometric Function Fit');
hold off

% 
% % Given that I have preferred vs null direction of neuron what to do next?
% % Choose range of ambiguous stimuli - based on fitted sigmoidal curve
% % Calculate choice eprobability of ambiguous stimuli