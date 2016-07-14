function play_flash(paramFile,stimFreq,stimDur,blockDur,tChar,minTR)

%% Displays a black/white full-field flicker
%
%   Usage:
%   play_flash(paramFile,stimFreq,stimDur,blockDur)
%
%   Inputs:
%   paramFile   - full path to output file containing timing of events, etc
%   stimFreq    - stimulus flicker frequency    (default = 16   [hertz])
%   stimDur     - duration of entire stimulus   (default = 336  [seconds])
%   blockDur    - duration of stimulus blocks   (default = 12   [seconds])
%   tChar           = {'t'}; % character(s) to signal a scanner trigger
%   minTR       - minimum time allowed between TRs (for use with recording triggers) (default = 0.25 [seconds])
%
%   Stimulus will flicker at 'stimFreq', occilating between flicker and
%   grey screen based on 'blockDur'
%
%   Written by Andrew S Bock Jul 2016

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
% stimulus frequency
if ~exist('stimFreq','var')
    stimFreq = 16; % seconds
end
% stimulus duration
if ~exist('stimDur','var')
    stimDur = 336; % seconds
end
% block duration
if ~exist('blockDur','var')
    blockDur = 12; % seconds
end
% scanner trigger
if ~exist('tChar','var') || isempty(tChar)
    tChar = {'t'};
end
% minimum time between TRs
if ~exist('minTR','var') || isempty(minTR)
    minTR = 0.25;
end
%% Save input variables
params.functionName     = mfilename;
params.gitInfo          = gitInfo;
params.userName         = userName;
params.stimFreq         = stimFreq;
params.stimDur          = stimDur;
params.blockDur         = blockDur;
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
%% Make images
greyScreen = grey*ones(fliplr(display.resolution));
blackScreen = black*ones(fliplr(display.resolution));
whiteScreen = white*ones(fliplr(display.resolution));
Texture(1) = Screen('MakeTexture', winPtr, blackScreen);
Texture(2) = Screen('MakeTexture', winPtr, whiteScreen);
Texture(3) = Screen('MakeTexture', winPtr, greyScreen);
%% Display Text, wait for Trigger
try
    commandwindow;
    Screen('FillRect',winPtr, grey);
    Screen('TextSize',winPtr,40);
    DrawFormattedText(winPtr, 'SCAN STARTING SOON, HOLD STILL!!!', ...
        'center',display.resolution(2)/3,[],[],[],[],[],0);
    Screen('Flip',winPtr);
    soundsc(sin(1:.5:1000)); % play 'ready' tone
    wait4T(keybs);  %wait for 't' from scanner.
    ListenChar(2);
    HideCursor;
    %% Drawing Loop
    breakIt = 0;
    frameCt = 0;
    TRct = 1;
    startTime = GetSecs;  %read the clock
    curFrame = 0;
    params.startDateTime    = datestr(now);
    params.endDateTime      = datestr(now); % this is updated below
    params.TRtime(TRct) = GetSecs;
    disp(['T ' num2str(TRct) ' received - 0 seconds']);
    lastT = startTime;
    elapsedTime = 0;
    while elapsedTime < stimDur && ~breakIt  %loop until 'esc' pressed or time runs out
        % get 't' from scanner
        [keyIsDown, secs, keyCode, ~] = KbCheck(-3);
        if keyIsDown % If *any* key is down
            % If 't' is one of the keys being pressed
            if sum(ismember(KbName(tChar),find(keyCode)))
                if (secs-lastT) > minTR
                    TRct = TRct + 1;
                    params.TRtime(TRct) = GetSecs;
                    disp(['T ' num2str(TRct) ' received - ' num2str(elapsedTime) ' seconds']);
                    lastT = secs;
                end
            end
        end
        % Flip between grey and flicker
        thisblock = floor(elapsedTime/blockDur);
        if mod(thisblock,2)
            if (elapsedTime - curFrame) > (1/(stimFreq*2))
                frameCt = frameCt + 1;
                Screen( 'DrawTexture', winPtr, Texture( mod(frameCt,2) + 1 )); % current frame
                % Flip to the screen
                Screen('Flip', winPtr);
                curFrame = GetSecs - startTime;
            end
        else
            Screen( 'DrawTexture', winPtr, Texture( 3 )); % current frame
            % Flip to the screen
            Screen('Flip', winPtr);
        end
        % update timers
        elapsedTime = GetSecs-startTime;
        params.endDateTime = datestr(now);
        % check to see if the "esc" button was pressed
        breakIt = escPressed(keybs);
    end
    sca;
    disp(['elapsedTime = ' num2str(elapsedTime)]);
    ListenChar(1);
    ShowCursor;
    Screen('CloseAll');
    %% Save params
    params.display = display;
    save(paramFile,'params');
catch ME
    Screen('CloseAll');
    ListenChar;
    ShowCursor;
    rethrow(ME);
end