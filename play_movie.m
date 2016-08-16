function play_movie(saveInfo,movieName,movieTime,indexisFrames,soundVol,tChar,TR,display)

%% Usage:
%
%   play_movie(paramFile,movieName,movieTime,indexisFrames,soundVol,tChar,TR,display)
%
%   Required inputs:
%   paramFile           - full path to output file containing timing of events, etc
%
%   Defaults:
%   movieName           - '/Users/abock/PixarShorts.mov'; % pixar shorts
%   movieTime           - [0 30]; % [start end] time of movie
%   endTime             - 30; % start time of movie
%   indexisFrames       - 0; % '1' changes the index from seconds to movie frames
%   soundVol            - 0; % can be a value between 0 (off) and 1 (full volume)
%   tChar               - {'t'}; % character(s) to signal a scanner trigger
%   TR                  - 0.8; % TR (seconds)
%
%   Written by Andrew S Bock Mar 2014

%% Set defaults
% Get git repository information
fCheck = which('GetGitInfo');
if ~isempty(fCheck)
    thePath = fileparts(mfilename('fullpath'));
    gitInfo = GetGITInfo(thePath);
else
    gitInfo = 'function ''GetGITInfo'' not found';
end
% Get user name
[~, userName] = system('whoami');
userName = strtrim(userName);
% movieName
if ~exist('movieName','var') || isempty(movieName)
    movieName = fullfile('/Users',userName,'PixarShorts.mov');
end
% movie start and end time
if ~exist('movieTime','var') || isempty(movieTime)
    movieTime = [0 30]; % [start end];
end
% time index (frames or seconds)
if ~exist('indexisFrames','var') || isempty(indexisFrames)
    indexisFrames = 0;
end
% sound volume
if ~exist('soundVol','var') || isempty(soundVol)
    soundVol = 0;
end
% scanner trigger
if ~exist('tChar','var') || isempty(tChar)
    tChar = {'t'};
end
% TR
if ~exist('TR','var') || isempty(TR)
    TR = 0.8;
end
% dispaly parameters
if ~exist('display','var') || isempty(display)
    display.distance = 106.5; % distance from screen (cm) - (UPenn - SC3T);
    display.width = 69.7347; % width of screen (cm) - (UPenn - SC3T);
    display.height = 39.2257; % height of screen (cm) - (UPenn - SC3T);
end
%% Save input variables
params.functionName     = mfilename;
params.gitInfo          = gitInfo;
params.userName         = userName;
params.subjectName      = saveInfo.subjectName;
params.TR               = TR;
params.scanDur          = movieTime(2) - movieTime(1);
params.movieName        = movieName;
params.movieTime        = movieTime;
params.indexisFrames    = indexisFrames;
params.soundVol         = soundVol;
%% Initial settings
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 2); % Skip sync tests
screens = Screen('Screens'); % get the number of screens
screenid = max(screens); % draw to the external screen
%% For Trigger
a = cd;
if a(1)=='/' % mac or linux
    a = PsychHID('Devices');
    for i = 1:length(a)
        d(i) = strcmp(a(i).usageName, 'Keyboard');
    end
    keybs = find(d);
else % windows
    keybs = [];
end
%% Define black and white
black = BlackIndex(screenid);
white = WhiteIndex(screenid);
grey = white/2;
%% Screen params
res = Screen('Resolution',max(Screen('screens')));
display.resolution = [res.width res.height];
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'UseRetinaResolution');
[winPtr] = PsychImaging('OpenWindow', screenid, grey);
[mint,~,~] = Screen('GetFlipInterval',winPtr,200);
display.frameRate = 1/mint; % 1/monitor flip interval = framerate (Hz)
display.screenAngle = pix2angle( display, display.resolution );
%% Play the movie
try
    commandwindow;
    % Display Text, wait for Trigger
    Screen('FillRect',winPtr, grey);
    Screen('TextSize',winPtr,40);
    DrawFormattedText(winPtr, 'SCAN STARTING SOON, HOLD STILL!!!', ...
        'center',display.resolution(2)/3,[],[],[],[],[],0);
    Screen('Flip',winPtr);
    ListenChar(2);
    HideCursor;
    % Setup movie
    [movieObj]=Screen('OpenMovie',winPtr,movieName);
    % Move to requested timeindex where texture loading should start:
    Screen('SetMovieTimeIndex', movieObj, movieTime(1), indexisFrames);
    Screen('PlayMovie',movieObj,1,0,soundVol);
    % Wait for trigger
    soundsc(sin(1:.5:1000)); % play 'ready' tone
    disp('Ready, waiting for trigger...');
    startTime = wait4T(tChar);  %wait for 't' from scanner.
    % Save params
    %TRct = 1; % first TRTime already recorded (first TR starts the movie)
    %disp(['T ' num2str(TRct) ' received - 0 seconds']);
    breakIt = 0;
    params.startDateTime    = datestr(now);
    params.endDateTime      = datestr(now); % this is updated below
    params.ScreenflipTime(1)  = GetSecs();
    %params.TRTime(1)    = GetSecs(); % first flipTime is PlayMovie
    lastFrame=-1;           % Presentation timestamp of last frame.
    frameCt=0;              % Number of loaded movie frames.
    flipct = 1; % first flipTime is PlayMovie
    disp(['Trigger received - ' params.startDateTime]);
    while GetSecs-startTime < (movieTime(2) - movieTime(1)) && ~breakIt  %loop until 'esc' pressed or time runs out
        % update timers
        elapsedTime = GetSecs-startTime;
        % check to see if the "esc" button was pressed
        breakIt = escPressed(keybs);
        %         % get 't' from scanner
        %         [keyIsDown, secs, keyCode, ~] = KbCheck(-3);
        %         if keyIsDown % If *any* key is down
        %             % If t is one of the keys being pressed
        %             if sum(ismember(KbName(tChar),find(keyCode)))
        %                 if (secs-params.TRTime(end)) > minTR
        %                     TRct = TRct+1;
        %                     params.TRTime(TRct) = secs;
        %                     disp(['T ' num2str(TRct) ' received - ' num2str(elapsedTime) ' seconds']);
        %                 end
        %             end
        %         end
        [frameTexture,frameTime] = Screen('GetMovieImage', winPtr, movieObj, 1);
        if (frameTexture>0 && frameTime>lastFrame)
            % Store its texture handle and exact movie timestamp in
            % arrays for later use:
            frameCt                         = frameCt + 1;
            params.frameTime(frameCt)       = frameTime;
            lastFrame=frameTime;
        else
            ShowCursor;
            ListenChar(1);
            return;
        end;
        Screen('DrawTexture',winPtr,frameTexture,[],[0 0 res.width res.height]);
        flipTime = Screen('Flip',winPtr);
        flipct = flipct+1;
        params.ScreenflipTime(flipct) = flipTime;
        Screen('Close',frameTexture);
        params.endDateTime = datestr(now);
        WaitSecs(0.001);
    end
    Screen('PlayMovie',movieObj,0); % Stop Movie
    Screen('CloseMovie',movieObj); % Close Move
    clear Screen;sca; % Clear the screen
    disp(['elapsedTime = ' num2str(elapsedTime)]);
    ListenChar(1);
    ShowCursor;
    Screen('CloseAll');
    %% Save params
    params.display = display;
    save(saveInfo.fileName,'params');
catch ME
    Screen('CloseAll');
    ListenChar(1);
    ShowCursor;
    rethrow(ME);
end