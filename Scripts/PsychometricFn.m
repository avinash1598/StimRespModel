% clc
% clear all
% close all

% expData = loadExpData('./Data/expData_g0_c1.mat');
% 
% trialMatrix1 = expData.trialMatrix;
% stimVector1 = trialMatrix1(:, 2);
% unqStimOrientations = unique(stimVector1)';
% 
% psychometricData = zeros(1, length(unqStimOrientations));
% 
% for i=1:length(unqStimOrientations)
%     stimOrientation = unqStimOrientations(i);
%     givenOrientationTrialIDxes = find(stimVector1 == stimOrientation);
%     decision = trialMatrix1(givenOrientationTrialIDxes, 4);
% 
%     percent_CCW = length(find(decision == -1)) / length(decision);
%     psychometricData(i) = percent_CCW;
% end
% 
% x = unqStimOrientations';
% y = psychometricData';
% 
% figure
% 
% % Sigmoid function
% sigmoidEq = fittype('1 / (1 + exp(-a*(x - b)))', 'independent', 'x', 'coefficients', {'a', 'b'});
% startPoint = [1, mean(x)];  % a = 1 (initial slope), b = mean(x) (midpoint)
% fitResult = fit(x, y, sigmoidEq, 'StartPoint', startPoint);
% 
% hold on
% scatter(x, y, 'b');
% plot(fitResult, 'b-');


expData = loadExpData('./Data/expData_g0.3_c1.mat');

trialMatrix1 = expData.trialMatrix;
stimVector1 = trialMatrix1(:, 2);
unqStimOrientations = unique(stimVector1)';

psychometricData = zeros(1, length(unqStimOrientations));

for i = 1:length(unqStimOrientations)
    stimOrientation = unqStimOrientations(i);
    givenOrientationTrialIDxes = find(stimVector1 == stimOrientation);
    decision = trialMatrix1(givenOrientationTrialIDxes, 4);

    percent_CCW = length(find(decision == -1)) / length(decision);
    psychometricData(i) = percent_CCW;
end

x = unqStimOrientations';
y = psychometricData';

% Sigmoid function
sigmoidEq = fittype('1 / (1 + exp(-a*(x - b)))', 'independent', 'x', 'coefficients', {'a', 'b'});
startPoint = [1, mean(x)];  % a = 1 (initial slope), b = mean(x) (midpoint)
fitResult = fit(x, y, sigmoidEq, 'StartPoint', startPoint);

figure
hold on
scatter(x, y, 'r');
plot(fitResult, 'r-');


expData = loadExpData('./Data/expData_g0.3_c0.1.mat');

trialMatrix1 = expData.trialMatrix;
stimVector1 = trialMatrix1(:, 2);
unqStimOrientations = unique(stimVector1)';

psychometricData = zeros(1, length(unqStimOrientations));

for i = 1:length(unqStimOrientations)
    stimOrientation = unqStimOrientations(i);
    givenOrientationTrialIDxes = find(stimVector1 == stimOrientation);
    decision = trialMatrix1(givenOrientationTrialIDxes, 4);

    percent_CCW = length(find(decision == -1)) / length(decision);
    psychometricData(i) = percent_CCW;
end

x = unqStimOrientations';
y = psychometricData';

% Sigmoid function
sigmoidEq = fittype('1 / (1 + exp(-a*(x - b)))', 'independent', 'x', 'coefficients', {'a', 'b'});
startPoint = [1, mean(x)];
fitResult = fit(x, y, sigmoidEq, 'StartPoint', startPoint);

scatter(x, y, 'g');
plot(fitResult, 'g-');

legend('', 'C1', '', 'C0.1');
xlabel("Orientation (deg)")
ylabel("% CCW")
title('Psychometric Function Fit');

hold off

% expData = loadExpData('./Data/expData_g0.3_c_unambi0.3.mat');
% 
% trialMatrix1 = expData.trialMatrix;
% stimVector1 = trialMatrix1(:, 2);
% unqStimOrientations = unique(stimVector1)';
% 
% psychometricData = zeros(1, length(unqStimOrientations));
% 
% for i = 1:length(unqStimOrientations)
%     stimOrientation = unqStimOrientations(i);
%     givenOrientationTrialIDxes = find(stimVector1 == stimOrientation);
%     decision = trialMatrix1(givenOrientationTrialIDxes, 4);
% 
%     percent_CCW = length(find(decision == -1)) / length(decision);
%     psychometricData(i) = percent_CCW;
% end
% 
% x = unqStimOrientations';
% y = psychometricData';
% 
% % Sigmoid function
% sigmoidEq = fittype('1 / (1 + exp(-a*(x - b)))', 'independent', 'x', 'coefficients', {'a', 'b'});
% startPoint = [1, mean(x)];
% fitResult = fit(x, y, sigmoidEq, 'StartPoint', startPoint);
% 
% scatter(x, y, 'k');
% plot(fitResult, 'k-');
% 
% legend('', 'G0-C1', '', 'G0.3-C1', '', 'G0.3-C_all_0.3', '', 'G0.3-C_unambi_0.3');
% xlabel("Orientation (deg)")
% ylabel("% CCW")
% title('Psychometric Function Fit');
% 
% hold off
