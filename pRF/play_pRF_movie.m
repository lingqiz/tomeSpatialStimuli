function play_pRF_movie(subj,runNum,images,TR,scanDur,display,redFrames)

%% Play pRF movie stimuli
%   'images' input created using 'make_bars.m'
%   Usage:
%   play_pRF_movie(subj,runNum,images,TR,scanDur,display,redFrames)
%
%   Written by Andrew S Bock Aug 2014

%% Set defaults
if ~exist('TR','var')
    TR = 0.8;
end
if ~exist('scanDur','var')
    scanDur = 336; % seconds
end
if ~exist('display','var')
    display.distance = 106.5; % distance from screen (cm) - (SC3T);  124.25 - (HUP6)
    display.width = 69.7347; % width of screen (cm) - (SC3T); 50.4 - (HUP6)
    display.skipChecks = 2;
    display.bkColor = [128 128 128];
    display.screenNum = max(Screen('Screens'));
end
% Make fixation dot color changes
if ~exist('redFrames','var')
    maxFrames = size(images,3);
    redFrames = zeros(1,maxFrames);
    minDiff = 0;
    while minDiff < ceil(4*TR*8); % dot color changes are separated by at least 4 TRs
        switches = sort(randperm(maxFrames,ceil(maxFrames/8/TR/20))); % ~every 20s
        minDiff = min(diff(switches));
    end
    ct = 0;
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
stim.redFrames = redFrames;
%% Initial settings
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 2); % Skip sync tests
screens = Screen('Screens'); % get the number of screens
screenid = max(screens); % draw to the external screen
%% For Trigger
a = cd;
if a(1)=='/' % mac or linux
    a = PsychHID('Devices');
    for i = 1:length(a), d(i) = strcmp(a(i).usageName, 'Keyboard'); end
    keybs = find(d);
else % windows
    keybs = [];
end
commandwindow
%% Define black and white
white = WhiteIndex(screenid);
black = BlackIndex(screenid);
grey = white/2;
%% Screen params
res = Screen('Resolution',max(Screen('screens')));
display.resolution = [res.width res.height];
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'UseRetinaResolution');
[winPtr, windowRect] = PsychImaging('OpenWindow', screenid, grey);
[mint,~,~] = Screen('GetFlipInterval',winPtr,200);
display.frameRate = 1/mint; % 1/monitor flip interval = framerate (Hz)
display.screenAngle = pix2angle( display, display.resolution );
%rect = Screen('Rect', winPtr );
[screenXpix, screenYpix] = Screen('WindowSize', winPtr);% Get the size of the on screen window
display.resolution = [screenXpix screenYpix];
[center(1), center(2)] = RectCenter(windowRect); % Get the center coordinate of the window
fix_mask = angle2pix(display,0.75); % For fixation mask (0.75 degree)
fix_dot = angle2pix(display,0.25); % For fixation cross (0.25 degree)
%% Dot stimulus params
% Set the blend function so that we get nice antialised edges
Screen('BlendFunction', winPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
%% stimulus specific params
for i = 1:size(images,3);
    tmp = images(:,:,i);
    Texture(i) = Screen('MakeTexture', winPtr, tmp);
end
try
    commandwindow;
    %% Display Text, wait for Trigger
    Screen('FillRect',winPtr, grey);
    Screen('TextSize',winPtr,40);
    DrawFormattedText(winPtr, 'SCAN STARTING SOON, HOLD STILL!!!', ...
        'center',display.resolution(2)/3,[],[],[],[],[],0);
    Screen('DrawDots', winPtr, [0;0], fix_dot,black, center, 1);
    Screen('Flip',winPtr);
    wait4T(keybs);  %wait for 't' from scanner.
    ListenChar(2);
    HideCursor;
    %% Drawing Loop
    breakIt = 0;
    Keyct = 0;
    Rct = 0;
    Gct = 0;
    curFrame = 1;
    curTR = 1;
    startTime = GetSecs;  %read the clock
    stim.startTime = startTime;
    stim.TRtime(curTR) = GetSecs;
    disp(['T ' num2str(curTR) ' received']);
    lastT = startTime;
    lastR = startTime;
    while GetSecs-startTime < scanDur && ~breakIt  %loop until 'esc' pressed or time runs out
        % update timers
        elapsedTime = GetSecs-startTime;
        % get 't' from scanner
        [keyIsDown, secs, keyCode, ~] = KbCheck(-3);
        if keyIsDown % If *any* key is down
            % If 't' is one of the keys being pressed
            if ismember(KbName('t'),find(keyCode))
                if (secs-lastT) > 0.25
                    curTR = curTR + 1;
                    stim.TRtime(curTR) = GetSecs;
                    disp(['T ' num2str(curTR) ' received']);
                    lastT = secs;
                end
            end
        end
        % Display 8 frames / TR
        if abs((elapsedTime / (TR / 8 )) - curFrame) > 0 %(TR / 6)
            curFrame = ceil( elapsedTime / (TR / 8 ));
        end
        % carrier
        Screen( 'DrawTexture', winPtr, Texture(curFrame)); % current frame
        % Fixation Mask
        Screen('FillOval',winPtr,grey,[screenXpix/2-fix_mask/2, ...
            screenYpix/2-fix_mask/2,screenXpix/2+fix_mask/2,screenYpix/2+fix_mask/2]);
        if redFrames(curFrame)
            Rct = Rct + 1;
            stim.RFrame(Rct) = GetSecs;
            Screen('DrawDots', winPtr, [0;0], fix_dot, [1 0 0], center, 1);
        else
            Gct = Gct + 1;
            stim.GFrame(Gct) = GetSecs;
            Screen('DrawDots', winPtr, [0;0], fix_dot, [0 1 0], center, 1);
        end
        % Flip to the screen
        Screen('Flip', winPtr);
        % record button presses
        [keyIsDown, secs, keyCode, ~] = KbCheck(-3);
        if keyIsDown % If *any* key is down
            % If r is one of the keys being pressed
            if ismember(KbName('r'),find(keyCode))
                if (secs-lastR) > 0.25
                    Keyct = Keyct+1;
                    lastR = secs;
                    stim.RT(Keyct) = GetSecs;
                    disp(['R ' num2str(Keyct) ' received']);
                end
            end
        end
        % check to see if the "esc" button was pressed
        breakIt = escPressed(keybs);
    end
    sca;
catch ME
    Screen('CloseAll');
    ListenChar;
    ShowCursor;
    rethrow(ME);
end
disp(['elapsedTime = ' num2str(elapsedTime)]);
ListenChar(1);
ShowCursor;
Screen('CloseAll');
if ~exist(fullfile('~/Desktop/MRI',[subj '_' date]),'dir')
    mkdir(fullfile('~/Desktop/MRI',[subj '_' date]));
end
save(fullfile('~/Desktop/MRI',[subj '_' date],['run' num2str(runNum)]),'stim');
return
