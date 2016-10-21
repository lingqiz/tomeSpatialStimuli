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
dbDir = ['/Users/' userName '/Dropbox-Aguirre-Brainard-Lab'];
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
tmp = datestr(now,2);
sessDate = tmp([1,2,4,5,7,8]);
%% Get the session name
sprintf(['\nSession Names:\n' ...
    '\n1 - session1_restAndStructure' ...
    '\n2 - session2_spatialStimuli' ...
    '\n3 - session3_OneLight\n'])
sessNum = input('Which session number? 1/2/3:\n');
if isempty(sessNum)
    error('no session number!');
end
%% Set output directory
if strcmp(userName,'connectome')
    outDir = fullfile(dbDir,'TOME_data',sessNames{sessNum},subjName,sessDate,'Stimuli');
else
    outDir = fullfile('/Users',userName,'TOME_data',sessNames{sessNum},subjName,sessDate,'Stimuli');
end
if ~exist(outDir,'dir')
    mkdir(outDir);
end
%% Get the stimulus type
sprintf(['\nRun Types:\n' ...
    '\n1 - MOVIE' ...
    '\n2 - RETINO' ...
    '\n3 - FLASH' ...
    '\n4 - FULL\n'])
runType = input('Which run type? 1/2/3/4:\n','s');
if isempty(runType)
    error('no run type!');
end
%% Get the run number
runNum = input('Run number?:\n','s');
if isempty(runNum)
    error('no run number!');
end
%% Get the function input
saveInfo.subjectName        = subjName;
switch runType
    case '1'
        movieName = fullfile(dbDir,'TOME_materials','stimulusFiles','PixarShorts.mov');
        switch runNum
            case '1'
                movieTime   = [1880 2216];
                runName     = 'tfMRI_MOVIE_AP_run01';
            case '2'
                movieTime   = [2216 2552];
                runName     = 'tfMRI_MOVIE_AP_run02';
            case '3'
                movieTime   = [892 1228];
                runName     = 'tfMRI_MOVIE_PA_run03';
            case '4'
                movieTime   = [1228 1564];
                runName     = 'tfMRI_MOVIE_PA_run04';
        end
        saveInfo.fileName   = fullfile(outDir,[runName '.mat']);
        % move file if re-running
        if exist(saveInfo.fileName,'file')
            if ~exist(fullfile(outDir,'abortedRuns'),'dir')
                mkdir(fullfile(outDir,'abortedRuns'));
            end
            system(['mv ' saveInfo.fileName ' ' fullfile(outDir,'abortedRuns',[runName '.mat'])]);
        end
        play_movie(saveInfo,movieName,movieTime);
    case '2'
        switch runNum
            case '1'
                runName     = 'tfMRI_RETINO_PA_run01';
            case '2'
                runName     = 'tfMRI_RETINO_PA_run02';
            case '3'
                runName     = 'tfMRI_RETINO_AP_run03';
            case '4'
                runName     = 'tfMRI_RETINO_AP_run04';
        end
        saveInfo.fileName   = fullfile(outDir,[runName '.mat']);
        % move file if re-running
        if exist(saveInfo.fileName,'file')
            if ~exist(fullfile(outDir,'abortedRuns'),'dir')
                mkdir(fullfile(outDir,'abortedRuns'));
            end
            system(['mv ' saveInfo.fileName ' ' fullfile(outDir,'abortedRuns',[runName '.mat'])]);
        end
        play_pRF(saveInfo);
    case '3'
        switch runNum
            case '1'
                runName     = 'tfMRI_FLASH_AP_run01';
            case '2'
                runName     = 'tfMRI_FLASH_PA_run02';
        end
        saveInfo.fileName   = fullfile(outDir,[runName '.mat']);
        % move file if re-running
        if exist(saveInfo.fileName,'file')
            if ~exist(fullfile(outDir,'abortedRuns'),'dir')
                mkdir(fullfile(outDir,'abortedRuns'));
            end
            system(['mv ' saveInfo.fileName ' ' fullfile(outDir,'abortedRuns',[runName '.mat'])]);
        end
        play_flash(saveInfo);
    case '4'
        movieName = fullfile(dbDir,'TOME_materials','stimulusFiles','WALL-E.mp4');
        switch runNum
            case '1'
                runName     = 'dMRI_T1w_T2w_run01';
            case '2'
                runName     = 'dMRI_T1w_T2w_run02';
        end
        movieStart = input('Movie start time (sec)? e.g. 0:\n');
        movieTime = [movieStart inf];
        saveInfo.fileName   = fullfile(outDir,[runName '.mat']);
        % move file if re-running
        if exist(saveInfo.fileName,'file')
            if ~exist(fullfile(outDir,'abortedRuns'),'dir')
                mkdir(fullfile(outDir,'abortedRuns'));
            end
            system(['mv ' saveInfo.fileName ' ' fullfile(outDir,'abortedRuns',[runName '.mat'])]);
        end
        play_movie(saveInfo,movieName,movieTime);
    otherwise
        disp('unknown run type');
end