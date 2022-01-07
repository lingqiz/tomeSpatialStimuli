%% Play the retinotopic stimulus

%% Add path
addpath(fullfile('.', 'pRF'));

%% Generate image stimuli
imFile = './data/imFile.mat'; % Output file that saves the images
imagesFull = make_bars(imFile); % Create the images that will be displayed

%% Save info
sID = 1; subName = 'LQZ';
fName = strcat(subName, '_', num2str(sID), '.mat');
savePath = fullfile('.', 'data', fName);

saveInfo.subjectName = subName;
saveInfo.fileName = savePath;

%% Play the stimulus
stimFile = load(fullfile('.', 'data', 'imFile.mat'));
play_pRF(saveInfo, stimFile.imagesFull);