function playStimulus

%   Master function that calls various stimulus fuctions:
%
%       play_movie
%       play_pRF
%       play_flash
%
%   Inputs to these functions are created through prompts to the user
%
%   Written by Andrew S Bock Jul 2016

%% Set defaults
% Get user name
[~, tmpName] = system('whoami');
userName = strtrim(tmpName);
% Set Dropbox directory
dbDir = ['/Users/' userName '/Dropbox (Aguirre-Brainard Lab)'];
disp(['Dropbox directory = ' dbDir]);
sessNames = {...
    'session1_restAndStructure' ...
    'session2_spatialStimuli' ...
    'session3_OneLight'};
%% Get the subject name
subjName = input('Subject name? e.g. TOME_3###:\n','s');
if isempty(subjName)
    error('no subject name!');
end
%% Get the session date
sessDate = input('Session date? e.g. 122516:\n','s');
if isempty(sessDate)
    error('no session date!');
end
%% Get the session name
sprintf(['\nSession Names:\n' ...
    '\n1 - session1_restAndStructure' ...
    '\n2 - session2_spatialStimuli' ...
    '\n3 - session3_OneLight\n'])
sessNum = input('Which session number? <1>/2/3:\n');
if isempty(sessNum)
    sessNum = 1;
end
%% Set output directory
outDir = fullfile(dbDir,'TOME_data',sessNames{sessNum},subjName,sessDate);
if ~exist(outDir,'dir')
    mkdir(outDir);
end
%% Get the stimulus type
runType = input('Which run type? movie/pRF/flash/DTI/T1/T2:\n','s');
if isempty(runType)
    error('no stimulus type!');
end
%% Get the run number
runNum = input('Run number?:\n','s');
if isempty(runNum)
    error('no run number!');
end
%% Get the function input
stimDir = fullfile(outDir,'StimulusFiles');
if ~exist(stimDir,'dir')
    mkdir(stimDir);
end
paramFile = fullfile(stimDir,[runType '_run' runNum '.mat']);
switch runType
    case 'movie'
        movieName = fullfile(dbDir,'TOME_materials','PixarShorts.mov');
        switch runNum
            case '1'
                movieTime = [1880 2216];
            case '2'
                movieTime = [2216 2552];
            case '3'
                movieTime = [892 1228];
            case '4'
                movieTime = [1228 1564];
        end
        play_movie(paramFile,movieName,movieTime);
    case 'pRF'
        play_pRF(paramFile,imagesFull,TR,scanDur,display,tChar,rChar,minTR)
    case 'flash'
        play_flash(paramFile,stimFreq,stimDur,blockDur,tChar,minTR);
    case DTI
        movieName = fullfile(dbDir,'TOME_materials','WALL-E.mov');
        switch runNum
            case '1'
                movieTime = [0 300];
            case '2'
                movieTime = [300 600];
            case '3'
                movieTime = [600 900];
            case '4'
                movieTime = [900 1200];
        end
        play_movie(paramFile,movieName,movieTime);
    case T1
        movieName = fullfile(dbDir,'TOME_materials','WALL-E.mov');
        movieTime = [0 300];
        play_movie(paramFile,movieName,movieTime);
    case T2
        movieName = fullfile(dbDir,'TOME_materials','WALL-E.mov');
        movieTime = [0 300];
        play_movie(paramFile,movieName,movieTime);
    otherwise
        disp('unknown run type');
end