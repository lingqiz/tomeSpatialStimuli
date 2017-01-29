function play_flash(saveInfo,stimFreq,scanDur,blockDur,tChar,TR,display)

%% Displays a black/white full-field flicker
%
%   Usage:
%   play_flash(paramFile,stimFreq,scanDur,blockDur,tChar,TR,display)
%
%   Required inputs:
%   paramFile           - full path to output file containing timing of events, etc
%
%   Defaults:
%   stimFreq            - stimulus flicker frequency    (default = 16   [hertz])
%   scanDur             - duration of entire stimulus   (default = 336  [seconds])
%   blockDur            - duration of stimulus blocks   (default = 12   [seconds])
%   tChar               - {'t'}; % character(s) to signal a scanner trigger
%   TR                  - 0.8; % TR (seconds)
%   display.distance    - 106.5; % distance from screen (cm) - (UPenn - SC3T);
%   display.width       - 69.7347; % width of screen (cm) - (UPenn - SC3T);
%   display.height      - 39.2257; % height of screen (cm) - (UPenn - SC3T);
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
    scanDur = 336; % seconds
end
% block duration
if ~exist('blockDur','var')
    blockDur = 12; % seconds
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
params.stimDur          = scanDur;
params.stimFreq         = stimFreq;
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
[winPtr, windowRect]            = PsychImaging('OpenWindow', screenid, grey);
[mint,~,~] = Screen('GetFlipInterval',winPtr,200);
display.frameRate = 1/mint; % 1/monitor flip interval = framerate (Hz)
display.screenAngle = pix2angle( display, display.resolution );
[center(1), center(2)]          = RectCenter(windowRect); % Get the center coordinate of the window
fix_dot                         = angle2pix(display,0.25); % For fixation cross (0.25 degree)
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
    Screen('DrawDots', winPtr, [0;0], fix_dot,black, center, 1);
    Screen('Flip',winPtr);
    ListenChar(2);
    HideCursor;
    soundsc(sin(1:.5:1000)); % play 'ready' tone
    disp('Ready, waiting for trigger...');
    startTime = wait4T(tChar);  %wait for 't' from scanner.
    %% Drawing Loop
    breakIt = 0;
    frameCt = 0;
    %TRct = 1;
    curFrame = 0;
    params.startDateTime    = datestr(now);
    params.endDateTime      = datestr(now); % this is updated below
    %params.TRtime(TRct) = GetSecs;
    %lastT = startTime;
    elapsedTime = 0;
    disp(['Trigger received - ' params.startDateTime]);
    while elapsedTime < scanDur && ~breakIt  %loop until 'esc' pressed or time runs out
        %         % get 't' from scanner
        %         [keyIsDown, secs, keyCode, ~] = KbCheck(-3);
        %         if keyIsDown % If *any* key is down
        %             % If 't' is one of the keys being pressed
        %             if sum(ismember(KbName(tChar),find(keyCode)))
        %                 if (secs-lastT) > minTR
        %                     TRct = TRct + 1;
        %                     params.TRtime(TRct) = GetSecs;
        %                     disp(['T ' num2str(TRct) ' received - ' num2str(elapsedTime) ' seconds']);
        %                     lastT = secs;
        %                 end
        %             end
        %         end
        % Flip between background and flicker
        thisblock = floor(elapsedTime/blockDur);
        if mod(thisblock,2)
            % flicker
            if (elapsedTime - curFrame) > (1/(stimFreq*2))
                frameCt = frameCt + 1;
                Screen( 'DrawTexture', winPtr, Texture( mod(frameCt,2) + 1 )); % current frame
                % Flip to the screen
                Screen('Flip', winPtr);
                curFrame = GetSecs - startTime;
            end
        else
            % background
            Screen( 'DrawTexture', winPtr, Texture( 1 )); % black screen
            % Flip to the screen
            Screen('Flip', winPtr);
        end
        % update timers
        elapsedTime = GetSecs-startTime;
        params.endDateTime = datestr(now);
        % check to see if the "esc" button was pressed
        breakIt = escPressed(keybs);
        WaitSecs(0.001);
    end
    sca;
    disp(['elapsedTime = ' num2str(elapsedTime)]);
    ListenChar(1);
    ShowCursor;
    Screen('CloseAll');
    %% Save params
    params.display = display;
    save(saveInfo.fileName,'params');
catch ME
    Screen('CloseAll');
    ListenChar;
    ShowCursor;
    rethrow(ME);
end