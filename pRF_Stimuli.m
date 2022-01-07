%% Play the retinotopic stimulus

% Add path
addpath(fullfile('.', 'pRF'));

% Define inputs
imFile = './data/imFile.mat'; % Output file that saves the images
imagesFull = make_bars(imFile); % Create the images that will be displayed

% Save info
savePath = './data/subData.mat';
saveInfo.subjectName = 'LQZ';
saveInfo.fileName = savePath;

% Play the stimulus
play_pRF(saveInfo, imagesFull);