function play_pRF(saveInfo, imagesFull, TR, scanDur, display, tChar, rChar)
%% Play pRF movie stimuli
%
%   Usage:
%   play_pRF(saveInfo,stimParams,TR,scanDur,display,tChar,rChar)
%
%   Required inputs:
%   saveInfo.subjectName    - subject name
%   saveInfo.fileName       - full path and name of output file
%
%   Defaults:
%   imagesFull              - image stimulus
%   TR                      - 0.8; TR (seconds)
%   scanDur                 - 336: scan duration (seconds)
%
%   UPenn - SC3T
%   display.distance        - 106.5; distance from screen (cm) 
%   display.width           - 69.7347; width of screen (cm) 
%   display.height          - 39.2257; height of screen (cm)
%
%   tChar                   - {'t'}; character(s) to signal a scanner trigger
%   rChar                   - {'r' 'g' 'b' 'y'}; character(s) to signal a button response
%
%   Written by Andrew S Bock Aug 2014

%% Set defaults parameters

% TR
if ~exist('TR','var') || isempty(TR)
    TR = 0.8;
end

% scan duration
if ~exist('scanDur','var') || isempty(scanDur)
    scanDur = 336; % seconds
end

% dispaly parameters
if ~exist('display','var') || isempty(display)
    display.distance = 106.5; % distance from screen (cm) - (UPenn - SC3T);
    display.width = 69.7347; % width of screen (cm) - (UPenn - SC3T);
    display.height = 39.2257; % height of screen (cm) - (UPenn - SC3T);
end

% scanner trigger
if ~exist('tChar','var') || isempty(tChar)
    tChar = {'t'};
end

% scanner trigger
if ~exist('rChar','var') || isempty(rChar)
    rChar = {'r' 'g' 'b' 'y'};
end

% Save input variables
params.functionName = mfilename;
params.subjectName = saveInfo.subjectName;
params.TR = TR;
params.scanDur = scanDur;

%% For Trigger
a = cd;
% mac or linux
if a(1) == '/'
    a = PsychHID('Devices');
    for i = 1:length(a)
        d(i) = strcmp(a(i).usageName, 'Keyboard');
    end
    keybs = find(d);
% windows
else
    keybs = [];
end

%% Make fixation dot color changes
if ~exist('redFrames','var') || isempty(redFrames)
    maxFrames = size(imagesFull,3);
    redFrames = zeros(1,maxFrames);
    minDiff = 0;
    maxDiff = inf;
    minTime = 5; % dot color changes are separated by at least 5s
    maxTime = 30; % dot color changes are separated no more than 30s

    while minDiff < (minTime * 8) || maxDiff > (maxTime * 8) % (8 frames / sec)
        switches = sort(randperm(maxFrames,ceil(maxFrames / (20 * 8) ))); % ~every 20s
        minDiff = min(diff(switches));
        maxDiff = max(diff(switches));
    end
    ct = 0;

    % Make a vector of 0's (green) and 1's (red), to use for chaging the
    %   color of the fixation dot
    for i = 1:length(switches)
        if ct
            if i ~= length(switches)
                redFrames(switches(i):switches(i+1)) = 1;
                ct = 0;
            else
                redFrames(switches(i):end) = 1;
                ct = 0;
            end
        else
            ct = ct + 1;
        end
    end
end
%% Initial settings
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 2); % Skip sync tests
screens = Screen('Screens'); % get the number of screens
screenid = max(screens); % draw to the external screen

%% Define black and white
white                           = WhiteIndex(screenid);
black                           = BlackIndex(screenid);
grey                            = white/2;

%% Screen params
res                             = Screen('Resolution',max(Screen('screens')));
display.resolution              = [res.width res.height];
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'UseRetinaResolution');
[winPtr, windowRect]            = PsychImaging('OpenWindow', screenid, grey);
[mint,~,~]                      = Screen('GetFlipInterval',winPtr,100);
display.frameRate               = 1/mint; % 1/monitor flip interval = framerate (Hz)
display.screenAngle             = pix2angle( display, display.resolution );
[screenXpix, screenYpix]        = Screen('WindowSize', winPtr);% Get the size of the on screen window
[center(1), center(2)]          = RectCenter(windowRect); % Get the center coordinate of the window
fix_mask                        = angle2pix(display,0.75); % For fixation mask (0.75 degree)
fix_dot                         = angle2pix(display,0.25); % For fixation cross (0.25 degree)

%% Dot stimulus params
% Set the blend function so that we get nice antialised edges
Screen('BlendFunction', winPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%% stimulus specific params
Texture = zeros(1, size(imagesFull, 3));
for i = 1:size(imagesFull, 3)
    tmp = imagesFull(:,:,i);
    Texture(i) = Screen('MakeTexture', winPtr, tmp);
end

% get rid of large variable to save space
clear tmp;
clear imagesFull;

%% Set to command window
commandwindow;

%% Run try/catch
try
    %% Display Text, wait for Trigger
    ListenChar(2);    
    soundsc(sin(1:.5:1000)); % play 'ready' tone
    disp('Ready, waiting for trigger...');
        
    % Wait for trigger
    target = 84; % 't' code
    code = 1;
    while ~ (length(code) == 1 && code == target)
        Screen('FillRect',winPtr, grey);
        Screen('DrawDots', winPtr, [0;0], fix_dot,black, center, 1);
        Screen('Flip',winPtr);

        [~, ~, keyCode, ~] = KbCheck;
        code = find(keyCode);
    end
    startTime = GetSecs;

    %% Drawing Loop
    breakIt                     = 0;
    Keyct                       = 0;
    Rct                         = 0;
    Gct                         = 0;
    curFrame                    = 1;
    params.startDateTime        = datestr(now);
    params.endDateTime          = datestr(now); % this is updated below
    params.GreenDotTime         = [];
    params.RedDotTime           = [];
    disp(['Trigger received - ' params.startDateTime]);
    while GetSecs-startTime < scanDur && ~breakIt  %loop until 'esc' pressed or time runs out
        % update timers
        elapsedTime = GetSecs-startTime;
        % check to see if the "esc" button was pressed
        breakIt = escPressed(keybs);
        % log button responses
        if CharAvail
            ch = GetChar;
            if sum(strcmp(ch,rChar))
                Keyct = Keyct+1;
                params.RT(Keyct) = GetSecs;
                disp(['Response ' num2str(Keyct) ' received']);
            end
            FlushEvents;
        end
        % Display 8 frames / TR
        if abs((elapsedTime / (TR / 8 )) - curFrame) > 0
            curFrame = ceil( elapsedTime / (TR / 8 ));
        end
        % carrier
        Screen( 'DrawTexture', winPtr, Texture(curFrame)); % current frame
        % Fixation Mask
        Screen('FillOval',winPtr,grey,[screenXpix/2-fix_mask/2, ...
            screenYpix/2-fix_mask/2,screenXpix/2+fix_mask/2,screenYpix/2+fix_mask/2]);
        if redFrames(curFrame)
            Rct = Rct + 1;
            params.RedDotTime(Rct) = GetSecs;
            Screen('DrawDots', winPtr, [0;0], fix_dot, [1 0 0], center, 1);
        else
            Gct = Gct + 1;
            params.GreenDotTime(Gct) = GetSecs;
            Screen('DrawDots', winPtr, [0;0], fix_dot, [0 1 0], center, 1);
        end
        % Flip to the screen
        Screen('Flip', winPtr);
        params.endDateTime      = datestr(now);
        WaitSecs(0.001);
    end
    sca;
    disp(['elapsedTime = ' num2str(elapsedTime)]);
    ListenChar(1);
    ShowCursor;
    Screen('CloseAll');

    %% Save params
    params.display = display;    
    save(saveInfo.fileName, 'params', '-v7.3');

catch ME
    Screen('CloseAll');
    ListenChar(1);
    ShowCursor;
    rethrow(ME);
end
