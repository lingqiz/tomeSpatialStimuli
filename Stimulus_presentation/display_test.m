function display_test(display)
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
red = [1 0 0];
green = [0 1 0];
blue = [0 0 1];
if ~exist('display','var')
    display.dist = 124.25; % distance from screen (cm) (HUP6)
    display.width = 50.4; % width of screen (cm) (HUP6)
    display.skipChecks = 2;
    display.bkColor = [128 128 128];
    display.screenNum = max(Screen('Screens'));
end
%display = OpenWindow(display);sca;
[window, windowRect] = PsychImaging('OpenWindow', screenid, grey);
        Screen('Flip', window);
try
    HideCursor;
    breakIt = 0;
    while ~breakIt
        [keyIsDown, secs, keyCode, ~] = KbCheck(-3);
        if keyIsDown % If *any* key is down
            % If t is one of the keys being pressed
            if (ismember(KbName('r'),find(keyCode)))
                [window, windowRect] = PsychImaging('OpenWindow', screenid, red);
                Screen('Flip', window);
            elseif (ismember(KbName('g'),find(keyCode)))
                [window, windowRect] = PsychImaging('OpenWindow', screenid, green);
                Screen('Flip', window);
            elseif (ismember(KbName('b'),find(keyCode)))
                [window, windowRect] = PsychImaging('OpenWindow', screenid, blue);
                Screen('Flip', window);
            elseif (ismember(KbName('k'),find(keyCode)))
                [window, windowRect] = PsychImaging('OpenWindow', screenid, black);
                Screen('Flip', window);
            elseif (ismember(KbName('w'),find(keyCode)))
                [window, windowRect] = PsychImaging('OpenWindow', screenid, white);
                Screen('Flip', window);
            elseif (ismember(KbName('y'),find(keyCode)))
                [window, windowRect] = PsychImaging('OpenWindow', screenid, grey);
                Screen('Flip', window);
            end
        end
        % check to see if the "esc" button was pressed
        breakIt = escPressed(keybs);
    end
    Screen('CloseAll');
    ListenChar(0);
    ShowCursor;
catch ME
    Screen('CloseAll');
    ListenChar(0);
    ShowCursor;
    rethrow(ME);
end