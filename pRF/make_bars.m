function [images,stimulus,params] = make_bars(params,imFile,paramsFile)

%% Makes drifting bar stimuli for pRF mapping.
%   Expected to be used in conjuction with play_pRF_movie
%
%   Usage:
%   [finalimages,images] = make_bars(params)
%
%   defaults:
%     params.screendims               = [69.85 39.37]; % cm  
%     params.resolution               = [1920 1200];
%     params.distance                 = 106.5; % cm 
%     params.step_nx                  = 16; % number of steps per one bar sweep
%     params.propScreen               = 1; % proportion of screen height (1 = full height of screen)
%     params.ringsize                 = 4; % proportion of screen radius (4 = 1/4 of radius)
%     params.motionSteps              = 8; % moving checks within bar
%     params.numSubRings              = 1;
%     params.display.backColorIndex   = 128; % set grey value
%     params.stimRgbRange             = [0 255];
%     params.numImages                = params.motionSteps*(params.step_nx); % total number of images in one pass (16 steps, 8 motion per step)
%     params.numReps                  = 3; % number of times to cycle through images
%     params.TR                       = 0.8; % TR in seconds
%
%   Written by Andrew S Bock Sep 2015

%% Set defaults (for 3T Prisma at UPenn)
if ~exist('params','var')
    params.screendims               = [69.85 39.37]; % cm  
    params.resolution               = [1920 1200];
    params.distance                 = 106.5; % cm 
    params.step_nx                  = 16; % number of steps per one bar sweep
    params.propScreen               = 1; % proportion of screen height (1 = full height of screen)
    params.ringsize                 = 4; % proportion of screen radius (4 = 1/4 of radius)
    params.motionSteps              = 8; % moving checks within bar
    params.numSubRings              = 1;
    params.display.backColorIndex   = 128; % set grey value
    params.stimRgbRange             = [0 255];
    params.numImages                = params.motionSteps*(params.step_nx); % total number of images in one pass (16 steps, 8 motion per step)
    params.numReps                  = 3; % number of times to cycle through images
    params.TR                       = 0.8; % TR in seconds
end
%% Pull out params
step_nx         = params.step_nx;                   % time for one pass (sec)
halfNumImages   = params.numImages./2;              % half of the images
numMotSteps     = params.motionSteps;               % 8, refers to the moving checks
numSubRings     = params.numSubRings;               % 1
bk              = params.display.backColorIndex;    % 128
minCmapVal      = min([params.stimRgbRange]);       % [0 255])
maxCmapVal      = max([params.stimRgbRange]);       % [0 255])
display.resolution  = params.resolution;            % [1920 1080]
display.distance    = params.distance;              % 106.5
display.width       = params.screendims(1);         % 69.85
display.height      = params.screendims(2);         % 69.85
outerRad            = 0.5*params.propScreen * (pix2angle(display,display.resolution(2),'height'));
ringWidth           = outerRad / params.ringsize;
if isfield(params, 'contrast')
    c = params.contrast;
    bg = (minCmapVal + maxCmapVal)/2;
    minCmapVal = round((1-c) * bg);
    maxCmapVal = round((1+c) * bg);
end
%% Initialize image template %%
disp('Creating images...');
m = round(angle2pix(display, 2*outerRad,'height'));
n = round(angle2pix(display, 2*outerRad,'width'));
[x,y]=meshgrid(linspace(-outerRad,outerRad,n),linspace(outerRad,-outerRad,m));
% here we crop the image if it is larger than the screen
% seems that you have to have a square matrix, bug either in my or
% psychtoolbox' code - so we make it square
if m>params.resolution(2)
    start  = round((m-params.resolution(2))/2);
    len    = params.resolution(2);
    y = y(start+1:start+len, start+1:start+len);
    x = x(start+1:start+len, start+1:start+len);
    m = len;
    n = len;
end;
% r = eccentricity;
r = sqrt (x.^2  + y.^2);
% loop over different orientations and make checkerboard
% first define which orientations
orientations    = deg2rad(0:45:360); % degrees -> rad
orientations    = orientations([1 6 3 8 5 2 7 4]); % shuffle order, remove 2pi
remake_xy       = zeros(1,params.numImages)-1;
remake_xy(1:length(remake_xy)/length(orientations):length(remake_xy)) = orientations;
original_x      = x;
original_y      = y;
% step size of the bar
step_x          = ((2*outerRad) - ringWidth/2) ./ step_nx;
step_startx     = -outerRad;%(step_nx-1)./2.*-step_x - (ringWidth./2);
% Loop that creates the one cycle images
tmpImages=zeros(m,n,halfNumImages*params.motionSteps,'uint8');
for imgNum=1:halfNumImages
    if remake_xy(imgNum) >=0
        x = original_x .* cos(remake_xy(imgNum)) - original_y .* sin(remake_xy(imgNum));
        y = original_x .* sin(remake_xy(imgNum)) + original_y .* cos(remake_xy(imgNum));
        % Calculate checkerboard.
        % Wedges alternating between -1 and 1 within stimulus window.
        % The computational contortions are to avoid sign=0 for sin zero-crossings
        wedges    = sign(round((cos((x+step_startx)*numSubRings*(2*pi/ringWidth)))./2+.5).*2-1);
        posWedges = find(wedges== 1);
        negWedges = find(wedges==-1);
        rings     = zeros(size(wedges));
        checks    = zeros(size(rings,1),size(rings,2),numMotSteps);
        for ii=1:numMotSteps,
            tmprings1 = sign(2*round((cos(y*numSubRings*(2*pi/ringWidth)+(ii-1)/numMotSteps*2*pi)+1)/2)-1);
            tmprings2 = sign(2*round((cos(y*numSubRings*(2*pi/ringWidth)-(ii-1)/numMotSteps*2*pi)+1)/2)-1);
            rings(posWedges) = tmprings1(posWedges);
            rings(negWedges) = tmprings2(negWedges);
            checks(:,:,ii)=minCmapVal+ceil((maxCmapVal-minCmapVal) * (wedges.*rings+1)./2);
        end;
        % reset starting point
        loX = step_startx - step_x;
    end;
    loX   = loX + step_x;
    hiX   = loX + ringWidth;
    % Ring completely disappear from view before it re-appears again in the middle.
    % we do this just be removing the second | from the window
    window = ( (x>=loX & x<=hiX) & r<outerRad);
    % yet another loop to be able to move the checks...
    tmpvar = zeros(m,n);
    tmpvar(window) = 1;
    tmpvar = repmat(tmpvar,[1 1 numMotSteps]);
    window = tmpvar == 1;
    img         = bk*ones(size(checks));
    img(window) = checks(window);
    tmpImages(:,:,(imgNum-1).*numMotSteps+1:imgNum.*numMotSteps) = uint8(img);
end
%% insert blank image
blankImage = uint8(ones(size(tmpImages,1),size(tmpImages,2)).*bk);
blankInd  = size(tmpImages,3)+1;
tmpImages(:,:,blankInd)   = blankImage;
%% Create final images
sweep = 1:params.numImages;
blank = ones(1,length(sweep)/2)*blankInd; % make a 'blank' bar sweep
% Make a single cycle through the images
tmpseq = [...
2*length(sweep) + sweep ...
blank ...
fliplr(1*length(sweep) + sweep) ...
fliplr(sweep) ...
blank ...
fliplr(3*length(sweep) + sweep) ...
fliplr(2*length(sweep) + sweep) ...
blank ...
1*length(sweep) + sweep ...
sweep ...
blank ...
3*length(sweep) + sweep ...
];
% Repeat the image cycle by 'params.numReps'
stimulus.seq = repmat(tmpseq,[1,params.numReps]);
images = tmpImages(:,:,stimulus.seq);
% stimulus timing
totalDur = (length(stimulus.seq) / params.motionSteps) * params.TR; % 8 frames / TR
stimulus.seqtiming = linspace(0,totalDur - (totalDur/length(stimulus.seq)),length(stimulus.seq));
disp('done.');
%% Save image and params files
save(imFile,'images');
svae(paramsFile,'params');