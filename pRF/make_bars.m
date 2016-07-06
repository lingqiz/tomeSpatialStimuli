function [imagesFull,images,stimulus,params] = make_bars(outFile,params)

%% Makes drifting bar stimuli for pRF mapping.
%   Expected to be used in conjuction with 'play_pRF'
%
%   Usage:
%   [imagesFull,images,stimulus,params] = make_bars(outFile,params)
%
%   Example:
%   outFile     = '/path/to/some/dir/outFile.mat'; % output image file
%   [images]    = make_bars(outFile);
%
%   defaults:
%     params.resolution               = [1920 1080]; % screen resolution
%     params.scanDur                  = 336; % scan duration in seconds
%     params.sweepDur                 = 16; % duration of a bar sweep in seconds
%     params.propScreen               = 1; % proportion of screen height (1 = full height of screen)
%     params.barsize                  = 1/4; % proportion of screen radius
%     params.checksize                = 2; % higher values = smaller checks
%     params.motionSteps              = 8; % moving checks within bar
%     params.display.backColorIndex   = 128; % set background grey value
%     params.stimRgbRange             = [0 255]; % color range
%     params.TR                       = 0.8; % TR in seconds
%     params.logBars                  = 0; % if = 1, bar width changes logarithmically with eccentricity
%     params.smallSteps               = 1; % if = 1, bar steps are smaller near the fovea
%
%   Written by Andrew S Bock Sep 2015

%% Set defaults (for 3T Prisma at UPenn)
if ~exist('params','var')
    params.resolution               = [1920 1080]; % screen resolution
    params.scanDur                  = 336; % scan duration in seconds
    params.sweepDur                 = 16; % duration of a bar sweep in seconds
    params.propScreen               = 1; % proportion of screen height (1 = full height of screen)
    params.barsize                  = 1/4; % proportion of screen radius
    params.checksize                = 2; % higher values = smaller checks
    params.motionSteps              = 8; % moving checks within bar
    params.display.backColorIndex   = 128; % set background grey value
    params.stimRgbRange             = [0 255]; % color range
    params.TR                       = 0.8; % TR in seconds
    params.logBars                  = 0; % if = 1, bar width changes logarithmically with eccentricity
    params.smallSteps               = 1; % if = 1, bar steps are smaller near the fovea
end
%% Pull out params
checkSize       = params.checksize;
sweepTRs        = params.sweepDur / params.TR;      % number of steps (TRs) for one pass
numMotSteps     = params.motionSteps;               % 8, refers to the moving checks
bk              = params.display.backColorIndex;    % 128
minCmapVal      = min([params.stimRgbRange]);       % [0 255])
maxCmapVal      = max([params.stimRgbRange]);       % [0 255])
outerRad        = 0.5*params.propScreen * (params.resolution(2));
if isfield(params, 'contrast')
    c = params.contrast;
    bg = (minCmapVal + maxCmapVal)/2;
    minCmapVal = round((1-c) * bg);
    maxCmapVal = round((1+c) * bg);
end
%% Calculate number of sweeps and buffer
numSweeps       = 8;                    % 4 orientations, 2 directions
numBlanks       = 4;                    % 4 blank periods of background
numImages       = numSweeps*(sweepTRs);  % bar orientations (8) x number of steps per bar pass
cycleDur        = numSweeps*params.sweepDur + numBlanks*params.sweepDur/2; % cycle duration (in seconds)
numCycles       = floor(params.scanDur / cycleDur);
bufferTime      = params.scanDur - numCycles*cycleDur;
%% Initialize image template %%
m = 2*outerRad;% height;
n = 2*outerRad; % width;
[original_x,original_y]=meshgrid(linspace(-outerRad,outerRad,n),linspace(outerRad,-outerRad,m));
% r = eccentricity;
r = sqrt (original_x.^2  + original_y.^2);
% define orientations
orientations        = deg2rad(0:45:360); % degrees -> rad
orientations        = orientations([1 6 3 8 5 2 7 4]); % shuffle order, remove 2pi
% Change x and y based on orientation (after each pass)
remake_xy           = -1*ones(1,numImages);
remake_xy(1:length(remake_xy)/length(orientations):length(remake_xy)) = orientations;
% step and ring size of the bar
if params.logBars
    % bar Width
    tmp = (logspace(0,1,sweepTRs/2+1)-1)/ max(logspace(0,1,sweepTRs/2+1)-1); % logspace between 0 and 1
    tmpWidths        = outerRad*diff(tmp);
    tmpWidths        = [fliplr(tmpWidths),tmpWidths];
    barWidths        = repmat(tmpWidths,1,numSweeps);
    % Window
    tmpL            = fliplr(outerRad - cumsum(tmpWidths));
    tmpR            = -(fliplr(tmpL));
else
    % check Width
    tmpWidths        = repmat(outerRad*params.barsize,1,sweepTRs);
    barWidths        = repmat(tmpWidths,1,numSweeps);
    % Window
    if params.smallSteps
        barWidth = outerRad*params.barsize;
        %         % Move bar using a quadratic
        %         tmp = ((1:sweepTRs/2-1).^2) / max((1:sweepTRs/2-1).^2);
        % Move bar logarithmically
        tmp = (logspace(0.1,1,sweepTRs/2-1)-1)/ max(logspace(0.1,1,sweepTRs/2-1)-1);
        % Steps based on leading edge
        LbarEdges = fliplr(-[0,tmp*(outerRad - barWidth)]);
        RbarEdges = [tmp*(outerRad - barWidth),outerRad];
        tmpR = [LbarEdges,RbarEdges]; % right edge
        tmpL = tmpR - barWidth; % left edge
    else
        tmpL            = (linspace(0,2*outerRad-outerRad*params.barsize,sweepTRs)) - outerRad;
        tmpR            = (linspace(outerRad*params.barsize,2*outerRad,sweepTRs)) - outerRad;
    end
end
lowX            = repmat(tmpL,1,numSweeps);
highX           = repmat(tmpR,1,numSweeps);
%% Create images
images=zeros(m,n,numImages*params.motionSteps,'uint8');
progBar = ProgressBar(numImages,'Creating images...');
for imgNum=1:numImages
    if remake_xy(imgNum) >=0
        % Rotate x and y based on bar orientation
        x = original_x .* cos(remake_xy(imgNum)) - original_y .* sin(remake_xy(imgNum));
        y = original_x .* sin(remake_xy(imgNum)) + original_y .* cos(remake_xy(imgNum));
    end
    % bars alternating between -1 and 1 within stimulus window.
    bars        = sign(round((cos( (checkSize*x-outerRad) * (2*pi/barWidths(imgNum)) ) ...
        )./2+.5).*2-1); % this last bit avoids bars == 0;
    posBars     = find(bars == 1);
    negBars     = find(bars ==-1);
    checks      = zeros(size(bars));
    allChecks   = zeros(size(checks,1),size(checks,2),numMotSteps);
    % Create checkerboard, move bars in opposite directions
    for ii=1:numMotSteps,
        posChecks = sign(2*round(( ...
            cos( checkSize*y*(2*pi/barWidths(imgNum) )...
            +(ii-1)/numMotSteps*2*pi)+1)/2)-1);
        negChecks = sign(2*round(( ...
            cos( checkSize*y*(2*pi/barWidths(imgNum) )...
            -(ii-1)/numMotSteps*2*pi)+1)/2)-1);
        checks(posBars) = posChecks(posBars);
        checks(negBars) = negChecks(negBars);
        allChecks(:,:,ii)=minCmapVal+ceil((maxCmapVal-minCmapVal) * (bars.*checks+1)./2);
    end
    % Make window to show checkboard
    tmpwindow = zeros(m,n);
    tmpwindow( (x>=lowX(imgNum) & x<=highX(imgNum)) & r<outerRad) = 1;
    window = logical(repmat(tmpwindow,[1 1 numMotSteps]));
    % Make images
    img         = bk*ones(size(allChecks));
    img(window) = allChecks(window);
    images(:,:,(imgNum-1)*numMotSteps+1:imgNum*numMotSteps) = uint8(img);
    progBar(imgNum);
end
%% insert blank image
blankImage = uint8(ones(size(images,1),size(images,2)).*bk);
blankInd  = size(images,3)+1;
images(:,:,blankInd)   = blankImage;
%% Create final sequence
sweep = 1:numImages;
bufferBlank = ones(1,params.motionSteps*(bufferTime / params.TR))*blankInd;
blank = ones(1,length(sweep)/2)*blankInd; % make a 'blank' bar sweep
% Make a single cycle through the images
tmpseq = [...
    blank ...
    1*length(sweep) + sweep ...         % diagonal (up-left)
    sweep ...                           % vertical (left-right)
    blank ...
    3*length(sweep) + sweep ...         % diagonal (up-right)
    2*length(sweep) + sweep ...         % horizontal (down)
    blank ...
    5*length(sweep) + sweep ...         % diagonal (down-right)
    4*length(sweep) + sweep ...         % vertical (right-left)
    blank ...
    7*length(sweep) + sweep ...         % diagonal (down-left)
    6*length(sweep) + sweep ...         % horizontal (up)
    ];
% Repeat the image cycle by 'params.numReps'
fullseq = repmat(tmpseq,[1,numCycles]);
stimulus.seq = [fullseq bufferBlank];
imagesFull = images(:,:,stimulus.seq);
% stimulus timing
totalDur = (length(stimulus.seq) / params.motionSteps) * params.TR; % 8 frames / TR
stimulus.seqtiming = linspace(0,totalDur - (totalDur/length(stimulus.seq)),length(stimulus.seq));
%% Save image and params files
disp('Saving output file...');
save(outFile,'images','stimulus','params','imagesFull','-v7.3');
disp('done.');