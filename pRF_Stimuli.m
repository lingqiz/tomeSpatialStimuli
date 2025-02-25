%% Play the retinotopic stimulus
% Add path
addpath(fullfile('.', 'pRF'));

%% Generate image stimuli (Run only once)
genStim = false;
if genStim
    imFile = './data/imFile.mat'; % Output file that saves the images
    imagesFull = make_bars(imFile); % Create the images that will be displayed
end

%% Save info
subName = 'PlaceHolder';
sID = input(sprintf('Subject %s, enter Session ID:', subName));

fName = strcat(subName, '_', num2str(sID), '.mat');
savePath = fullfile('.', 'data', fName);

saveInfo.subjectName = subName;
saveInfo.fileName = savePath;

% Play the stimulus
stimFile = load(fullfile('.', 'data', 'imFile.mat'));
stimFile = stimFile.imagesFull;

play_pRF(saveInfo, stimFile);