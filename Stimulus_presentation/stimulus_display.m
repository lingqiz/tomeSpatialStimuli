function stimulus_display(type,subj,run,blockdur,display)
%
%   stimulus_display_fresh(subj,run,display,type)
%
%   Written by Andrew S Bock Aug 2014

%% Subject and Scan info
if ~exist('type','var')
    type = 'dots';
end
if ~exist('subj','var')
    subj = 'foo';
end
if ~exist('run','var')
    run = 999;% if *test* only 2 TRs
end
if ~exist('blockdur','var')
    blockdur = 16; % duration of blocks, in seconds
end
% stimulus type
if strcmp(type,'dots')
    disps = {'Grey' 'Dots_left' 'Dots_right'};
elseif strcmp(type,'checker')
    disps = {'Grey' 'HemiChecker_left' 'HemiChecker_right'};
elseif strcmp(type,'MP')
    disps = {'Grey' 'Fullfield_M' 'Fullfield_P'}; %
else
    error('type not recognized');
end
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
if ~exist('display','var')
    display.dist = 124.25; % distance from screen (cm) (HUP6)
    display.width = 50.4; % width of screen (cm) (HUP6)
    display.skipChecks = 2;
    display.bkColor = [128 128 128];
    display.screenNum = max(Screen('Screens'));
end
%display = OpenWindow(display);sca;
res = Screen('Resolution',max(Screen('screens')));
display.resolution = [res.width res.height];
[winPtr, windowRect] = PsychImaging('OpenWindow', screenid, grey);
[mint,~,~] = Screen('GetFlipInterval',winPtr,200);
display.frameRate = 1/mint; % 1/monitor flip interval = framerate (Hz)
display.screenAngle = pix2angle( display, display.resolution );

rect = Screen('Rect', winPtr );
[screenXpix, screenYpix] = Screen('WindowSize', winPtr);% Get the size of the on screen window
display.resolution = [screenXpix screenYpix];
[center(1), center(2)] = RectCenter(windowRect); % Get the center coordinate of the window
%% Run params
% De Bruijn sequence
% Make sure Grey (1) always comes first, and is also the last in the sequence.
DeBruijn_sequence(1).vals = [1,2,1,1,3,3,2,2,3,1,2,1,1,3,3,2,2,3,1];
DeBruijn_sequence(2).vals = [1,3,3,2,1,1,2,2,3,1,3,3,2,1,1,2,2,3,1];
DeBruijn_sequence(3).vals = [1,2,2,3,1,3,3,2,1,1,2,2,3,1,3,3,2,1,1];
DeBruijn_sequence(4).vals = [1,3,2,3,3,1,2,2,1,1,3,2,3,3,1,2,2,1,1];
DeBruijn_sequence(5).vals = [1,3,2,2,3,3,1,1,2,1,3,2,2,3,3,1,1,2,1];
DeBruijn_sequence(6).vals = [1,3,1,1,2,2,3,3,2,1,3,1,1,2,2,3,3,2,1];
DeBruijn_sequence(7).vals = [1,2,2,3,2,1,3,3,1,1,2,2,3,2,1,3,3,1,1];
DeBruijn_sequence(8).vals = [1,2,3,3,2,2,1,3,1,1,2,3,3,2,2,1,3,1,1];
DeBruijn_sequence(9).vals = [1,3,3,2,1,2,2,3,1,1,3,3,2,1,2,2,3,1,1];
DeBruijn_sequence(10).vals = [1,2,2,1,3,2,3,3,1,1,2,2,1,3,2,3,3,1,1];
DeBruijn_sequence(999).vals = [1,2,1,1,3,3,2,2,3,1];
if run == 999
    DeBruijn = DeBruijn_sequence(run).vals;
elseif run > 10
    DeBruijn = DeBruijn_sequence(mod(run,10)).vals;
else
    DeBruijn = DeBruijn_sequence(run).vals;
end
%% Dot stimulus params
dpix = (angle2pix(display,12))/display.frameRate; % number of pixels to move for each frame (12 degree/s)
dotSizePixels = angle2pix(display,0.5); % Set the dot size in pixels (0.1 degree)
fix_dot = angle2pix(display,0.25); % For fixation cross (0.25 degree)
fix_mask = angle2pix(display,0.75); % For fixation mask (0.75 degree)
% Outward
outnumDots = 300;
outxStart = (rand(1, outnumDots) -.5) .* screenXpix;
outyStart = (rand(1, outnumDots) -.5) .* screenYpix;
% Inward
innumDots = 300;
inxStart = (rand(1, innumDots) -.5) .* screenXpix;
inyStart = (rand(1, innumDots) -.5) .* screenYpix;
edgeindx = [ones(1,screenYpix-10)*-screenXpix/2+10,(-screenXpix/2+10):0, ...
    1:screenXpix/2-10,ones(1,screenYpix-10)*screenXpix/2-10,fliplr(0:screenXpix/2-10), ...
    fliplr(-screenXpix/2+10:-1)];
edgeindy = [-screenYpix/2+10:0,1:screenYpix/2-10,ones(1,screenXpix-10)*screenYpix/2-10, ...
    fliplr(0:screenYpix/2-10),fliplr(-screenYpix/2+10:-1),ones(1,screenXpix-10)*-screenYpix/2];
% Random
randnumDots = 200;
randxStart = (rand(1, randnumDots) -.5) .* screenXpix;
randyStart = (rand(1, randnumDots) -.5) .* screenYpix;
randslope = rand(1,randnumDots)*(2*pi);
% Alternate between black-in, white-out, and black-out, white-in
colorct = 0;
% Set the blend function so that we get nice antialised edges
Screen('BlendFunction', winPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
if strcmp(type,'checker') || strcmp(type,'MP')
    for B = 1:length(disps)
        stim{B}.StimulusType = disps{B};
        if ~strcmp('Grey',stim{B}.StimulusType)
            %% default stimulus params
            % may be overwritten for specific experiments
            % carrier
            stim{B}.cyclePerDeg = .5; % checkerboard
            stim{B}.flickerHz = 8;
            % stim
            stim{B}.radDeg = 6; % circular aperature
            stim{B}.widthDeg = 2; % bar width
            stim{B}.randAng = 1; % logical - if 0 bar orientation varies in steps of 45 deg, if 1 orient is random
            stim{B}.degPerSecRange = [ 1 2 ]; % range of speeds
            if ~stim{B}.randAng, stim{B}.degPerSecRange = [1 1]*stim{B}.widthDeg; end % non overlapping steps
            stim{B}.durBlanksSec = 15; % durantion of blanks between bar sweeps
            stim{B}.nSweepsBetweenBlanks = 4; % number of sweeps separating two blank periods
            stim{B}.begendBlankSec = 5; % duration of blanks at beginning and end of run
            % time
            stim{B}.secPerTR = 2;
            stim{B}.numTRs = 8;
            stim{B}.durSec = stim{B}.secPerTR*stim{B}.numTRs; %
            % others
            stim{B}.fixMaskRadDeg = 0.25; % mask around fixation (0.25 is taken up by the simon display)
            %stim.seed = 8194; % always the same 'random' bar orient (Zach: 8194 best seed out of 1:1e5)
            %% stimulus specific params
            switch stim{B}.StimulusType
                case 'Fullfield_P'
                    stim{B}.flickerHz = 0.5;
                    % stimulus aperture
                    stim{B}.blockdur = 1;
                    stim{B}.Fs = 30;
                    dt=1/stim{B}.Fs;t=(0:dt:stim{B}.blockdur-dt)';
                    stim{B}.onecycle = [ones(stim{B}.durSec/stim{B}.flickerHz,1);-ones(stim{B}.durSec/stim{B}.flickerHz,1)];
                    stim{B}.fulltc = repmat(stim{B}.onecycle,stim{B}.durSec/stim{B}.flickerHz,1);
                    stim{B}.centersDeg = zeros(size(stim{B}.fulltc,1),2);
                    stim{B}.centersDeg(stim{B}.fulltc == 0) = 999; % so that blanks -> out of the screen == not shown
                    stim{B}.radDeg = min(display.screenAngle)/2;
                    stim{B}.squares = 'ecc';
                case 'Fullfield_M'
                    stim{B}.flickerHz = 16;
                    stim{B}.cutfreq = 0.5; % cycles/degree cuffoff for use with blur.m
                    % stimulus aperture
                    stim{B}.blockdur = 1;
                    stim{B}.Fs = 30;
                    dt=1/stim{B}.Fs;t=(0:dt:stim{B}.blockdur-dt)';
                    stim{B}.onecycle = cos(2*pi*stim{B}.blockdur*t);
                    stim{B}.fulltc = repmat(stim{B}.onecycle,stim{B}.durSec/stim{B}.flickerHz,1);
                    stim{B}.centersDeg = zeros(size(stim{B}.fulltc,1),2);
                    stim{B}.centersDeg(stim{B}.fulltc == 0) = 999; % so that blanks -> out of the screen == not shown
                    stim{B}.radDeg = min(display.screenAngle)/2;
                    stim{B}.squares = 'ecc';
                case {'HemiChecker','HemiChecker_left','HemiChecker_right'}
                    stim{B}.flickerHz = 8;
                    stim{B}.cutfreq = 0.5; % cycles/degree cuffoff for use with blur.m
                    % stimulus aperture
                    stim{B}.blockdur = 1;
                    stim{B}.Fs = 30;
                    dt=1/stim{B}.Fs;t=(0:dt:stim{B}.blockdur-dt)';
                    stim{B}.onecycle = cos(2*pi*stim{B}.blockdur*t);
                    stim{B}.fulltc = repmat(stim{B}.onecycle,stim{B}.durSec/stim{B}.flickerHz,1);
                    stim{B}.centersDeg = zeros(size(stim{B}.fulltc,1),2);
                    stim{B}.centersDeg(stim{B}.fulltc == 0) = 999; % so that blanks -> out of the screen == not shown
                    stim{B}.radDeg = min(display.screenAngle)/2;
                    stim{B}.squares = 'square';
            end
            stim{B}.nFrames = ceil( stim{B}.durSec * display.frameRate );
            if any( 2*stim{B}.radDeg > display.screenAngle )
                error( 'resize stimulus -- you are outside the bounds of your monitor' )
                Screen('CloseAll'); % ListenChar(1);
            end
            % carrier
            ndrfitpos = 10;
            stim_tmp = stim; stim_tmp.widthDeg = 100;%max(display.screenAngle)/2;
            if strcmp(stim{B}.StimulusType,'Fullfield_P')
                ph = 1;
                % 360 deg wheel
                %carrier = makeFullfieldImage(display, stim{B});
                carrier = checker(0,display,16,32,.241);
                img1 = carrier;
                img2 = -carrier+1;
                %                 img1 = 255*carrier;
                %                 img2 = 255*(-carrier+1);
                %             img1 = rescale2RGB( carrier, 1 ); % rescale img to RGB vals (100% contrast)
                %             img2 = rescale2RGB( -carrier, 1 ); % phase-reversed
                
                tmp1 = zeros(size(img1,1),size(img1,2),3);
                tmp2 = tmp1;
                
                tmp1(:,:,1) = img1;
                tmp1(:,:,2) = img2;
                tmp2(:,:,2) = img1;
                tmp2(:,:,1) = img2;
                
                stim{B}.img(1).img = tmp1;
                stim{B}.img(2).img = tmp2;
                
                Texture{B}(1,ph) = Screen('MakeTexture', winPtr, tmp1);
                Texture{B}(2,ph) = Screen('MakeTexture', winPtr, tmp2);
            elseif strcmp(stim{B}.StimulusType,'Fullfield_M')
                ph = 1;
                % 360 deg wheel
                stim{B}.ang = 9;
                stim{B}.rad = 2;
                carrier = checker(0,display,6,12,.165);
                carrier(carrier == 0) = -1; % rescale2RGB needs range -1 to 1
                carrier = 0.25*carrier; % make more grey
                carrier = rescale2RGB(carrier,1);
                carrier = blur(0,display,carrier,stim{B}.cutfreq); %.
                diff = 128 - mean(carrier(:));
                carrier = carrier + diff; % center at 128
                for m = 1:stim{B}.Fs
                    stim{B}.img(m).img = round(128 + (128 - carrier)*stim{B}.onecycle(m));
                    stim{B}.img(m).img(stim{B}.img(m).img>255) = 255;
                    Texture{B}(m,ph) = Screen('MakeTexture', winPtr, stim{B}.img(m).img);
                end
            elseif strcmp(stim{B}.StimulusType,'HemiChecker') || ...
                    strcmp(stim{B}.StimulusType,'HemiChecker_left') || ...
                    strcmp(stim{B}.StimulusType,'HemiChecker_right')
                ph = 1;
                % 360 deg wheel
                stim{B}.ang = 9;
                stim{B}.rad = 2;
                carrier = checker(0,display,12,24,.218);
                carrier(carrier == 0) = -1; % rescale2RGB needs range -1 to 1
                %carrier = 0.25*carrier; % make more grey
                carrier = rescale2RGB(carrier,1);
                %carrier = blur(0,display,carrier,stim{B}.cutfreq); %.
                diff = 128 - mean(carrier(:));
                carrier = carrier + diff; % center at 128
                for m = 1:stim{B}.Fs
                    stim{B}.img(m).img = round(128 + (128 - carrier)*stim{B}.onecycle(m));
                    stim{B}.img(m).img(stim{B}.img(m).img>255) = 255;
                    Texture{B}(m,ph) = Screen('MakeTexture', winPtr, stim{B}.img(m).img);
                end
            end
            clear stim_tmp
            Mask = zeros(display.resolution);
            % reformat stim{B}.centersDeg in rect space
            stimCentersPix = angle2pix( display, stim{B}.centersDeg );
            destRect = repmat( rect, size(stim{B}.fulltc,1), 1 ) + [stimCentersPix stimCentersPix ];
        end
    end
end
try
    commandwindow;
    %% Display Text, wait for Trigger
    Screen('TextSize',winPtr,40);
    DrawFormattedText(winPtr, 'SCAN STARTING SOON, HOLD STILL!!!', ...
        'center',display.resolution(2)/3,[],[],[],[],[],0);
    Screen('DrawDots', winPtr, [0;0], fix_dot, black, center, 1);
    Screen('Flip',winPtr);
    wait4T(keybs);  %wait for 't' from scanner.
    ListenChar(2);
    HideCursor;
    %% Drawing Loop
    for B = 1:length(DeBruijn)
        blockname = disps(DeBruijn(B));
        breakIt = 0;
        frameCnt = 1;
        timeStamp = 0;
        startTime = GetSecs;  %read the clock
        % Blinking fixation dot
        blink(1) = rand(1)*(blockdur/3);
        blink(2) = (rand(1)*(blockdur/3))+(blockdur/3);
        blink(3) = (rand(1)*(blockdur/3))+2*(blockdur/3);
        % Vary between flashing fixation, and flashing in hemifield for dots
        dtct = [round(rand(1)),round(rand(1)),round(rand(1))];
        % Random location for hemifield red dots
        if strcmp(blockname,'Dots_left')
            dtlc = [randi(screenXpix/2),randi(screenYpix); ...
                randi(screenXpix/2),randi(screenYpix); ...
                randi(screenXpix/2),randi(screenYpix)];
        else
            dtlc = [randi(screenXpix/2)+screenXpix/2,randi(screenYpix); ...
                randi(screenXpix/2)+screenXpix/2,randi(screenYpix); ...
                randi(screenXpix/2)+screenXpix/2,randi(screenYpix)];
        end
        stim{DeBruijn(B)}.keys(B).blink = blink;
        stim{DeBruijn(B)}.keys(B).RTime = 0;
        Keyct = 0;
        blinkDate = [];
        if strcmp(blockname,'Grey')
            startTime = GetSecs;
            while ~breakIt && GetSecs - startTime < blockdur
                % Fill screen with grey
                Screen('FillRect',winPtr,grey,...
                    [0,0,screenXpix,screenYpix]);
                % Fixation Cross
                Screen('DrawDots', winPtr, [0;0], fix_dot, black, center, 1);
                % Blinking dot
                if floor(blink(1)) < GetSecs-startTime && GetSecs-startTime < ceil(blink(1))-.85;
                    blinkDate(1) = now;
                    if strcmp(type,'dots')
                        Screen('DrawDots', winPtr, [0;0], fix_dot, [1 0 0], center, 1);
                    else
                        Screen('DrawDots', winPtr, [0;0], fix_dot, white, center, 1);
                    end
                elseif floor(blink(2)) < GetSecs-startTime && GetSecs-startTime < ceil(blink(2))-.85;
                    blinkDate(2) = now;
                    if strcmp(type,'dots')
                        Screen('DrawDots', winPtr, [0;0], fix_dot, [1 0 0], center, 1);
                    else
                        Screen('DrawDots', winPtr, [0;0], fix_dot, white, center, 1);
                    end
                elseif floor(blink(3)) < GetSecs-startTime && GetSecs-startTime < ceil(blink(3))-.85;
                    blinkDate(3) = now;
                    if strcmp(type,'dots')
                        Screen('DrawDots', winPtr, [0;0], fix_dot, [1 0 0], center, 1);
                    else
                        Screen('DrawDots', winPtr, [0;0], fix_dot, white, center, 1);
                    end
                end
                
                % Flip to the screen
                Screen('Flip', winPtr);
                breakIt = escPressed(keybs);
            end
        elseif strcmp(blockname,'Dots_left') || strcmp(blockname,'Dots_right')
            startTime = GetSecs;
            colorct = colorct +1;
            framect = 0;
            drawred = 0;
            redindx = randi(randnumDots);
            while ~breakIt && GetSecs - startTime < blockdur  % Stimulus drawing loop (exits when any button is pressed)
                % Draw the dots
                if ~mod(colorct,2)
                    Screen('DrawDots', winPtr, [outxStart; outyStart], dotSizePixels, white, center, 1);
                    Screen('DrawDots', winPtr, [inxStart; inyStart], dotSizePixels, black, center, 1);
                else
                    Screen('DrawDots', winPtr, [outxStart; outyStart], dotSizePixels, black, center, 1);
                    Screen('DrawDots', winPtr, [inxStart; inyStart], dotSizePixels, white, center, 1);
                end
                % Random Dots
                Screen('DrawDots', winPtr, [randxStart(1:randnumDots/2); randyStart(1:randnumDots/2)], dotSizePixels, white, center, 1);
                Screen('DrawDots', winPtr, [randxStart(randnumDots/2:randnumDots); randyStart(randnumDots/2:randnumDots)], dotSizePixels, black, center, 1);
                % Make one dot red (for subject attention/detection task)
                if drawred
                    Screen('DrawDots', winPtr, [xpos(outnumDots+innumDots+redindx);...
                        ypos(outnumDots+innumDots+redindx)], dotSizePixels, [1 0 0], center, 1);
                end
                % Hemifield mask
                if strcmp(blockname,'Dots_left')
                    Screen('FillRect',winPtr,grey,...
                        [screenXpix/2,0,screenXpix,screenYpix]);
                elseif strcmp(blockname,'Dots_right')
                    Screen('FillRect',winPtr,grey,...
                        [0,0,screenXpix/2,screenYpix]);
                end
                % Fixation Mask
                Screen('FillOval',winPtr,grey,[screenXpix/2-fix_mask/2, ...
                    screenYpix/2-fix_mask/2,screenXpix/2+fix_mask/2,screenYpix/2+fix_mask/2]);
                % Fixation Cross
                Screen('DrawDots', winPtr, [0;0], fix_dot, black, center, 1);
                % Flash red dot either on fixation cross, or in visual
                % hemifield
                if floor(blink(1)) < GetSecs-startTime && GetSecs-startTime < ceil(blink(1))-.85;
                    blinkDate(1) = now;
                    if dtct(1)
                        Screen('DrawDots', winPtr, [0;0], dotSizePixels, [1 0 0], ...
                            [dtlc(1,1),dtlc(1,2)], 1);
                    else
                        Screen('DrawDots', winPtr, [0;0], fix_dot, [1 0 0], center, 1);
                    end
                elseif floor(blink(2)) < GetSecs-startTime && GetSecs-startTime < ceil(blink(2))-.85;
                    blinkDate(2) = now;
                    if dtct(2)
                        Screen('DrawDots', winPtr, [0;0], dotSizePixels, [1 0 0], ...
                            [dtlc(2,1),dtlc(2,2)], 1);
                    else
                        Screen('DrawDots', winPtr, [0;0], fix_dot, [1 0 0], center, 1);
                    end
                elseif floor(blink(3)) < GetSecs-startTime && GetSecs-startTime < ceil(blink(3))-.85;
                    blinkDate(3) = now;
                    if dtct(3)
                        Screen('DrawDots', winPtr, [0;0], dotSizePixels, [1 0 0], ...
                            [dtlc(3,1),dtlc(3,2)], 1);
                    else
                        Screen('DrawDots', winPtr, [0;0], fix_dot, [1 0 0], center, 1);
                    end
                end
                % Flip to the screen
                Screen('Flip', winPtr);
                framect = framect + 1;
                % Find the line connecting each dot to the center of the screen, move the
                % dots along this line using dpix.
                %% Outward
                outslope = outyStart./outxStart;
                outtmppos = outxStart > 0; % right side of screen
                outtmpneg = outxStart < 0; % left side of screen
                outnewX(outtmppos) = outxStart(outtmppos) + sqrt((dpix^2)./(1+outslope(outtmppos).^2));
                outnewX(outtmpneg) = outxStart(outtmpneg) - sqrt((dpix^2)./(1+outslope(outtmpneg).^2));
                outxStart = outnewX;
                outyStart = outnewX.*outslope;
                % Find the dots that are off the screen, move them back to the
                % center (or at least very close, as setting them to [0,0] breaks
                % the code)
                outind = outxStart>screenXpix/2 | outyStart>screenYpix/2 ...
                    | outxStart<-screenXpix/2 | outyStart<-screenYpix/2;
                outxStart(outind) = (rand(1,sum(outind)) -.5) .* screenXpix/100;
                outyStart(outind) = (rand(1,sum(outind)) -.5) .* screenYpix/100;
                %% Inward
                inslope = inyStart./inxStart;
                intmppos = inxStart > 0; % right side of screen
                intmpneg = inxStart < 0; % left side of screen
                innewX(intmppos) = inxStart(intmppos) - sqrt((dpix^2)./(1+inslope(intmppos).^2));
                innewX(intmpneg) = inxStart(intmpneg) + sqrt((dpix^2)./(1+inslope(intmpneg).^2));
                inxStart = innewX;
                inyStart = innewX.*inslope;
                % Find the dots at center, move them to a random position
                inind = sqrt(inxStart.^2 + inyStart.^2)<fix_mask/2;
                %     inxStart(inind) = (rand(1,sum(inind)) -.5) .* screenXpix;
                %     inyStart(inind) = (rand(1,sum(inind)) -.5) .* screenYpix;
                tmpind = find(inind);
                for tmpi = 1:length(tmpind)
                    newind = randi(length(edgeindx));
                    inxStart(tmpind(tmpi)) = edgeindx(newind);
                    inyStart(tmpind(tmpi)) = edgeindy(newind);
                end
                %% Random
                %randslope = randyStart./randxStart;
                randxStart = randxStart + (dpix*sin(randslope));
                randyStart = randyStart - (dpix*cos(randslope));
                % Find the dots that are off the screen, move them back to the
                % center (or at least very close, as setting them to [0,0] breaks
                % the code)
                randind = randxStart>screenXpix/2 | randyStart>screenYpix/2 ...
                    | randxStart<-screenXpix/2 | randyStart<-screenYpix/2;
                randxStart(randind) = (rand(1,sum(randind)) -.5) .* screenXpix;
                randyStart(randind) = (rand(1,sum(randind)) -.5) .* screenYpix;
                %% Reset the "drawred" variable
                %     if ismember(redindx,find(randind))
                %         drawred = 0;
                %         redindx = randi(randnumDots);
                %     end
                %     if ~drawred && ~mod(framect,180);
                %         drawred = 1;
                %     end
                %%
                breakIt = escPressed(keybs);
                
            end
        elseif strcmp(blockname,'Fullfield_P') || strcmp(blockname,'Fullfield_M') || ...
                strcmp(blockname,'HemiChecker') || ...
                strcmp(blockname,'HemiChecker_left') || ...
                strcmp(blockname,'HemiChecker_right')
            imgblock = stim{DeBruijn(B)};
            texblock = Texture{DeBruijn(B)};
            breakIt = 0;
            frameCnt = 1;
            timeStamp = 0;
            display.timeStamp = zeros( 1, imgblock.nFrames );
            startTime = GetSecs;  %read the clock
            blink(1) = rand(1)*(blockdur/3);
            blink(2) = (rand(1)*(blockdur/3))+(blockdur/3);
            blink(3) = (rand(1)*(blockdur/3))+2*(blockdur/3);
            stim{DeBruijn(B)}.keys(B).blink = blink;
            stim{DeBruijn(B)}.keys(B).RTime = 0;
            Keyct = 0;
            blinkDate = [];
            while GetSecs-startTime < blockdur && ~breakIt  %loop until 'esc' pressed or time runs out
                % update timers
                elapsedTime = GetSecs-startTime;
                curTR = ceil( elapsedTime / imgblock.secPerTR );
                if strcmp(blockname,'Fullfield_M') || ...
                        strcmp(blockname,'HemiChecker') || ...
                        strcmp(blockname,'HemiChecker_left') || ...
                        strcmp(blockname,'HemiChecker_right')
                    curHz = ceil( elapsedTime * imgblock.Fs );
                    if curHz > imgblock.Fs
                        curHz = imgblock.Fs;
                    end
                end
                flickState =  round(mod( elapsedTime * imgblock.flickerHz, 1)) +1;
                if strcmp(blockname,'Fullfield_M') || ...
                        strcmp(blockname,'HemiChecker') || ...
                        strcmp(blockname,'HemiChecker_left') || ...
                        strcmp(blockname,'HemiChecker_right')
                    flickState =  round(mod(elapsedTime * imgblock.flickerHz,1)*imgblock.Fs)+1;
                    if flickState > imgblock.Fs
                        flickState = imgblock.Fs;
                    end
                end
                driftState = 1;
                if strcmp(blockname,'Fullfield_P')
                    % carrier
                    Screen( 'DrawTexture', winPtr, texblock(flickState,driftState), [], destRect(curTR,:)); % curTR
                elseif strcmp(blockname,'Fullfield_M') || ...
                        strcmp(blockname,'HemiChecker') || ...
                        strcmp(blockname,'HemiChecker_left') || ...
                        strcmp(blockname,'HemiChecker_right')
                    driftState = 1;
                    % carrier
                    Screen( 'DrawTexture', winPtr, texblock(flickState,driftState), [], destRect(curHz,:)); % cutHz
                    % Hemifield mask
                    if strcmp(blockname,'HemiChecker_left')
                        Screen('FillRect',winPtr,grey,...
                            [screenXpix/2,0,screenXpix,screenYpix]);
                    elseif strcmp(blockname,'HemiChecker_right')
                        Screen('FillRect',winPtr,grey,...
                            [0,0,screenXpix/2,screenYpix]);
                    end
                end
                % Fixation Mask
                Screen('FillOval',winPtr,grey,[screenXpix/2-fix_mask/2, ...
                    screenYpix/2-fix_mask/2,screenXpix/2+fix_mask/2,screenYpix/2+fix_mask/2]);
                % Fixation Cross
                Screen('DrawDots', winPtr, [0;0], fix_dot, black, center, 1);
                if floor(blink(1)) < GetSecs-startTime && GetSecs-startTime < ceil(blink(1))-.85;
                    blinkDate(1) = now;
                    Screen('DrawDots', winPtr, [0;0], fix_dot, white, center, 1);
                elseif floor(blink(2)) < GetSecs-startTime && GetSecs-startTime < ceil(blink(2))-.85;
                    blinkDate(2) = now;
                    Screen('DrawDots', winPtr, [0;0], fix_dot, white, center, 1);
                elseif floor(blink(3)) < GetSecs-startTime && GetSecs-startTime < ceil(blink(3))-.85;
                    blinkDate(3) = now;
                    Screen('DrawDots', winPtr, [0;0], fix_dot, white, center, 1);
                end
                % Flip to the screen
                Screen('Flip', winPtr);
                display.timeStamp(frameCnt) = GetSecs-startTime;
                % record button presses
                [keyIsDown, secs, keyCode, ~] = KbCheck(-3);
                if keyIsDown % If *any* key is down
                    % If t is one of the keys being pressed
                    if (ismember(KbName('r'),find(keyCode)))
                        if (secs-stim{DeBruijn(B)}.keys(B).RTime(end)) > 0.25
                            Keyct = Keyct+1;
                            stim{DeBruijn(B)}.keys(B).RTime(Keyct) = secs;
                            stim{DeBruijn(B)}.keys(B).RDate(Keyct) = now;
                            fprintf('*** R received ***\n');
                        end
                    end
                end
                stim{DeBruijn(B)}.keys(B).blinkDate = blinkDate;
                % check to see if the "esc" button was pressed
                breakIt = escPressed(keybs);
                frameCnt = frameCnt + 1;
            end
        end
    end
    sca;
catch ME
    disp(num2str(B));
    Screen('CloseAll');
    ListenChar;
    ShowCursor;
    rethrow(ME);
end
ListenChar(1);
ShowCursor;
Screen('CloseAll');
if ~exist(fullfile('~/Desktop/MRI',[subj '_' date]),'dir')
    mkdir(fullfile('~/Desktop/MRI',[subj '_' date]));
end
save(fullfile('~/Desktop/MRI',[subj '_' date],['run' num2str(run)]),'stim','-v7.3');
return
