% clc
% clear all
% close all

% TODO: set direction of neurons to none

rng('shuffle');

expData = loadExpData('expData_corr.mat');

trialMatrix = expData.trialMatrix;
trialResponses = expData.trialResponses;
neuronPrefOrientations = expData.preferredOrientation;

ambiguousStimuliOrientation = 90; % In degrees
neuronPreferredDirection = (neuronPrefOrientations > 90) * (-1) + (neuronPrefOrientations <= 90) * (1); % What about neurons which have 90 degree orientation
ambiguousStimTrialsIDx = find(trialMatrix(:, 2) == ambiguousStimuliOrientation);
ambiguousStimTrials = (trialMatrix(:, 2) == ambiguousStimuliOrientation);

nNeurons = length(neuronPrefOrientations);
neuronsCPs = zeros(nNeurons, 2); % Stores choice probability and its significance

ambiguousTrialsResponses = trialResponses(ambiguousStimTrials, :, :);

% Find noise responses for each neuron
trlSpkCnt = sum(ambiguousTrialsResponses, 3);
meanSpkCntAcrossTrials = mean(trlSpkCnt, 1);
sdDevSpkCntAcrossTrials = std(trlSpkCnt, 1);

zScoredTrlSpkCnt = (trlSpkCnt - meanSpkCntAcrossTrials)./sdDevSpkCntAcrossTrials;
corrMatrix = zeros(nNeurons, nNeurons);

for i=1:nNeurons
    for j=1:nNeurons
        if i~=j
            [r, p] = corr(zScoredTrlSpkCnt(:,i), zScoredTrlSpkCnt(:,j));
            corrMatrix(i, j) = r; % Take into account the significance of corr val as well
        end
    end
end

[sortedOrientations, sortIdx] = sort(neuronPrefOrientations, 'ascend');
sortedCorrMatrix = corrMatrix(sortIdx, sortIdx);

%%

figure
sgt = sgtitle('With noise correlation');

% subplot(2, 1, 1)
imagesc(sortedCorrMatrix);
idx = (1:15:nNeurons);
xticks(idx);
yticks(idx);
xticklabels(round(sortedOrientations(idx)))
yticklabels(round(sortedOrientations(idx)))
cb = colorbar;
cb.Label.String = 'Noise correlation coeff';
axis square;
xlabel("Preferred orientations (deg)")
ylabel("Preferred orientations (deg)")
% title("Noise correlation matrix")
set(gca, 'YDir', 'normal');  % Set the origin to the bottom-left corner