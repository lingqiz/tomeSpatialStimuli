% Stimuli MASTER script
%
%   This script provides details for running various functions that
%   display visual stimuli, specifically:
%
%   Retinotopic stimulus
%   Movie stimulus
%   Light flux stimulus
%
%   For each of the functions below, use 'help <functionName>' for
%   information to modify the default parameters
%
%   Written by Andrew S Bock Jul 2016

%% Retinotopic stimulus
% Define inputs
imFile          = '/data/project/subject/session/<someName>.mat'; % Output file that saves the images
imagesFull      = make_bars(imFile); % Create the images that will be displayed
paramFile       = '/data/project/subject/session/<someName>.mat'; % Output file that saves various parameters

% Play the stimulus
play_pRF(paramFile,imagesFull);
%% Movie stimulus
% Define inputs
movieName       = '/path/to/PixarShorts.mov'; % Define the movie to be shown
paramFile       = '/data/project/subject/session/<someName>.mat'; % Output file that saves various parameters
movieTime       = [0 336]; % Define the [start stop] time for the movie (seconds)

% Play the stimulus
play_movie(paramFile,movieName,movieTime);
%% Light flux stimulus
% Define inputs
paramFile       = '/data/project/subject/session/<someName>.mat'; % Output file that saves various parameters
stimFreq        = 16;   % flicker frequency (seconds)
stimDur         = 336;  % total stimulus duration (seconds)
blockDur        = 12;   % duration of stimulus block (seconds)

% Play the stimulus
play_flash(paramFile);