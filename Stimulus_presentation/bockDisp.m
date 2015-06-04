function bockDisp


%% Initial settings
AssertOpenGL;
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
minCmapVal = black; 
maxCmapVal = white;
%% Screen params
if ~exist('display','var')
    display.dist = 180.2; % distance from screen (cm) (7T)
    display.width = 69.7347; % width of screen (cm) (7T)
    display.height = 39.2257;
    display.skipChecks = 2;
    display.bkColor = grey;
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
% to allow blending
Screen('BlendFunction', winPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
%% Required params
outerRad = pix2angle(display,min(display.resolution))/2; % radius of visual display area
barWidth = outerRad/4; % 1/4 of visual display area - used for bar width
cycle = 600;
numMotSteps = 8; % checkerboard motion within the bar
TR = 2;
tempFreq = 2;
stimframe = 1./tempFreq./numMotSteps;
stimFrames = cycle./(1./tempFreq./numMotSteps);
blankDur = 14; % seconds of blank screen
initBlank = 4; % initial blank period, useful to make timing correct (i.e. exactly 10 minutes)
%% Initialize image template %%%
m = round(2*angle2pix(display,outerRad));
n = round(2*angle2pix(display,outerRad));
[x,y]=meshgrid(linspace(-outerRad,outerRad,n),linspace(outerRad,-outerRad,m));
% here we crop the image if it is larger than the screen
    % seems that you have to have a square matrix, bug either in my or
    % psychtoolbox' code - so we make it square
    if m>display.resolution(2),
        start  = round((m-display.resolution(2))/2);
        len    = display.resolution(2);
        y = y(start+1:start+len, start+1:start+len);
        x = x(start+1:start+len, start+1:start+len);
        m = len;
        n = len;
    end;
% eccentricity
r = sqrt (x.^2  + y.^2);
% orientations
orientations = deg2rad((0:45:360));
orientations = orientations([1 6 3 8 5 2 7 4]); % definitely need to remove the last (9th) value, not sure why the shuffle
%orientations = orientations([1 5 2 6 3 7 4 8]); % definitely need to remove the last (9th) value, not sure why the shuffle
numImages = cycle*numMotSteps/length(orientations); % should an even number
halfNumImages = numImages./2;
remake_xy = zeros(1,numImages) - 1;
remake_xy(1:length(remake_xy)/length(orientations):length(remake_xy)) = orientations;
% step size of the bar
step_nx      = cycle./TR/8;
step_x       = (2*outerRad) ./ step_nx;
step_startx  = (step_nx-1)./2.*-step_x - (barWidth./2);
% not really used
softmask = ones(m);
% Create Images
%images=zeros(m,n,halfNumImages*motionSteps,'uint8');
images=zeros(m,n,halfNumImages*numMotSteps);
original_x   = x;
original_y   = y;
progBar = ProgressBar(halfNumImages,['Creating ' num2str(halfNumImages) ' images']);
for imgNum=1:halfNumImages
    % Create a checkerboard, oriented according to the value found in
    % remake_xy
    if remake_xy(imgNum) >=0
        x = original_x .* cos(remake_xy(imgNum)) - original_y .* sin(remake_xy(imgNum));
        y = original_x .* sin(remake_xy(imgNum)) + original_y .* cos(remake_xy(imgNum));
        % Calculate checkerboard.
        % Wedges alternating between -1 and 1 within stimulus window.
        % The computational contortions are to avoid sign=0 for sin zero-crossings
        wedges    = sign(round((cos((x+step_startx)*(2*pi/barWidth)))./2+.5).*2-1);
        posWedges = find(wedges== 1);
        negWedges = find(wedges==-1);
        rings     = zeros(size(wedges));
        checks    = zeros(size(rings,1),size(rings,2),numMotSteps);
        for ii=1:numMotSteps,
            tmprings1 = sign(2*round((cos(y*(2*pi/barWidth)+(ii-1)/numMotSteps*2*pi)+1)/2)-1);
            tmprings2 = sign(2*round((cos(y*(2*pi/barWidth)-(ii-1)/numMotSteps*2*pi)+1)/2)-1);
            rings(posWedges) = tmprings1(posWedges);
            rings(negWedges) = tmprings2(negWedges);
            checks(:,:,ii)=minCmapVal+ceil((maxCmapVal-minCmapVal) * (wedges.*rings+1)./2);
        end;        
        % reset starting point
        loX = step_startx - step_x;
    end;
    % Define area of screen to present bar (i.e. window)
    loX   = loX + step_x;
    hiX   = loX + barWidth;
    window = ( (x>=loX & x<=hiX) & r<outerRad);
    tmpvar = zeros(m,n);
    tmpvar(window) = 1;
    % Copy this window for all motion steps (e.g. 8 steps)
    tmpvar = repmat(tmpvar,[1 1 numMotSteps]);
    window = tmpvar == 1;
    % Create background for all motion steps (e.g. 8 steps)
    img         = grey*ones(size(checks));
    % Create checkerboard bar in the bar window
    img(window) = checks(window);
    %images(:,:,(imgNum-1).*numMotSteps+1:imgNum.*numMotSteps) = uint8(img); %#ok<*BDSCA>
    % Copy img into images matrix
    images(:,:,(imgNum-1)*numMotSteps+1:imgNum*numMotSteps) = img;
    progBar(imgNum);
end
%% pad outside of square image to fit full screen
edges = (max(display.resolution) - size(images,1))/2;
tmp = ones(display.resolution(2),display.resolution(1),size(images,3))*grey;
tmp(:,edges+1:edges+size(images,1),:) = images;
images = tmp;
images = rescale2RGB(images,1);
clear tmp;
%% insert blanks
blankImage= ones(m,n,1)*grey;
blankInd  = size(images,3)+1;
images(:,:,blankInd)   = blankImage;
%% Show the stimulus
[winPtr, windowRect] = PsychImaging('OpenWindow', screenid, grey);
for i=1:size(images,3);
    Texture(i) = Screen('MakeTexture',winPtr,images(:,:,i));
end

startTime = GetSecs;  %read the clock
txCt = 1;
tic
while GetSecs-startTime<32;
    elapsedTime = GetSecs-startTime;
    curTR = ceil(elapsedTime/imgblock.secPerTR*numMotSteps);
    Screen( 'DrawTexture', winPtr, Texture(curTR)); % cutHz
    % Flip to the screen
    Screen('Flip', winPtr);
end
toc