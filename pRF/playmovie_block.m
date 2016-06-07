function playmovie_block(block,fix,subj,run,fromTime,nBlocks,blockDur,soundVol,indexisFrames,moviename)

%% Usage:
%
%   playmovie_block([,block],[,fix],[,subject],[,run],[,fromTime],[,nBlocks],
%   [,blockDur],[,soundVol],[,indexisFrames],[,moviename])
%
%   defaults:
%       block = 0; % just move (block = 1; movie and gray screen blocks)
%       fix = 0; % no fixation cross (fix = 1; fixation cross present)
%       subject = 'test' % subject name
%       run = 999 % runNum (e.g. 1)
%       fromTime = 0; % start time of movie
%       nBlocks = 1; % number of times to loop
%       blockDur = 60; % duration of block (in seconds)
%       soundVol = 0; % can be a value between 0 (off) and 1 (full volume)
%       indexisFrames = 0;
%       moviename = '/Users/abock/Downloads/Raiders/Raiders.mp4';
%
%   Note: currently coded to record triggers for TRs > 0.5s.
%
%   Written by Andrew S Bock Mar 2014
%% Initial settings
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 2); % Skip sync tests
%% keyboard variables
%   Switch KbName into unified mode: It will use the names of the OS-X
%   platform on all platforms in order to make this script portable:
KbName('UnifyKeyNames');
esc=KbName('ESCAPE');
%% Stimulus variables
if nargin<1
    block = 0;
end
if nargin<2
    fix = 0;
end
if nargin<3
    subj = 'test';
end
if nargin<4
    run = 999;
end
if nargin<5
    fromTime = 0;
end
if nargin<6
    nBlocks = 1;
end
if nargin<7
    blockDur = 60;
end
if nargin<8
    soundVol = 0;
end
if nargin<9
    indexisFrames = 0; % If set to 1 then the fromTime and toTime parameters are
    % interpreted as frameindex (starting with 0 for the first frame), instead
    % of seconds.
end
if nargin<10
    %moviename = '/Users/abock/Desktop/MRI/The_Artist_2011_1080p/The_Artist_2011_1080p.mp4';
    %moviename = '/Users/abock/Downloads/Raiders/Raiders.mp4';
    moviename = '/Users/abock/Downloads/ALL.mov';
    %moviename = [PsychtoolboxRoot 'PsychDemos/MovieDemos/DualDiscs.mov'];
    %moviename = '/Users/abock/Downloads/The_Artist_2011_1080p/The_Artist_2011_1080p.mp4';
    %moviename = '/Users/abock/Desktop/The_Artist_2011_1080p/The_Artist_2011_1080p_flip_lr.mov';
end

if block == 0
    blocktype = 'just movie';
else
    blocktype = 'movie and gray screen';
end
if fix == 1
    fixtype = 'fixation cross present';
else
    fixtype = 'no fixation cross present';
end
data.moviename = moviename;
data.fromTime = fromTime;
data.nBlocks = nBlocks;
data.blockDur = blockDur;
data.soundVol = soundVol;
data.blocktype = blocktype;
data.fixtype = fixtype;
%% Set up Initial Screen

res = Screen('Resolution',max(Screen('screens')));

s.size = [10,20,24,5,28]; % pixels sizes
s.gap = 3; %pixels
s.c = ceil([res.width res.height]/2);
for i=1:length(s.size)
    s.rect{i} = [-s.size(i)/2+s.c(1),-s.size(i)/2+s.c(2),+s.size(i)/2+s.c(1),+s.size(i)/2+s.c(2)];
end

[win,~]=Screen('OpenWindow',max(Screen('screens')),0);
background = [128,128,128]; % Set background color
Screen('FillRect',win, background);
Screen('TextSize',win,40);
DrawFormattedText(win, 'SCAN STARTING SOON, HOLD STILL!!!', ...
    'center',res.height/3,[],[],[],[],[],0);
Screen('Flip',win);
HideCursor;
commandwindow;
%ListenChar(2);
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
%%
[movie]=Screen('OpenMovie',win,moviename);
wait4T(keybs);  %wait for 't' from scanner.
TRct = 1; % first TRTime already recorded (first TR starts the movie)
data.startTime = GetSecs();
data.startDate = now;
data.flipTime(1) = data.startTime;
data.flipDate(1) = data.startDate;% first flipTime is PlayMovie
data.TRTime(1) = data.startTime; % first flipTime is PlayMovie
data.TRDate(1) = data.startDate;
try
    for B = 1:nBlocks
        tic
        %% Play movie
        if B == 1
            data.blockTimes = now;
        else
            data.blockTimes = [data.blockTimes now];
        end
        % Move to requested timeindex where texture loading should start:
        Screen('SetMovieTimeIndex', movie, fromTime+blockDur/2*(B-1), indexisFrames);
        lastpts=-1;          % Presentation timestamp of last frame.
        count=0;            % Number of loaded movie frames.
        Screen('PlayMovie',movie,1,0,soundVol);
        if B == 1
            flipct = 1; % first flipTime is PlayMovie
        end
        if block == 0;
            while (toc < blockDur)
                [keyIsDown, secs, keyCode, ~] = KbCheck(-3);
                if keyIsDown % If *any* key is down
                    % If t is one of the keys being pressed
                    if (ismember(KbName('t'),find(keyCode)))
                        if (secs-data.TRTime(end)) > 0.5
                            TRct = TRct+1;
                            data.TRTime(TRct) = secs;
                            data.TRDate(TRct) = now;
                            fprintf('*** T received ***\n');
                        end
                    elseif (ismember(esc,find(keyCode)))
                        Screen('PlayMovie',movie,0); % Stop Movie
                        Screen('CloseMovie',movie); % Close Move
                        clear Screen; % Clear the screen
                        Screen('CloseAll');sca; % Close the screen
                        ShowCursor;
                        ListenChar(1);
                        return
                    end
                end
                [movietexture,pts] = Screen('GetMovieImage', win, movie, 1);
                if (movietexture>0 && pts>lastpts)
                    % Store its texture handle and exact movie timestamp in
                    % arrays for later use:
                    count=count + 1;
                    data.texids(count)=movietexture;
                    data.texpts(count)=pts;
                    lastpts=pts;
                else
                    ShowCursor;
                    ListenChar(1);
                    return;
                end;
                Screen('DrawTexture',win,movietexture,[],[0 0 res.width res.height]);
                if fix == 1
                    Screen('FillOval',win,[0,0,0],s.rect{5}); % background of outer fixation circle
                    Screen('FillOval',win,[255,255,255],s.rect{3}); % background of inner fixation circle
                    Screen('FillRect',win,[0,0,0],[s.c(1)-s.size(2)/2,s.c(2)-s.gap/2,s.c(1)+s.size(2)/2,s.c(2)+s.gap/2]); %horizontal bar in cross
                    Screen('FillRect',win,[0,0,0],[s.c(1)-s.gap/2,s.c(2)-s.size(2)/2,s.c(1)+s.gap/2,s.c(2)+s.size(2)/2]); %vertical bar in cross
                    Screen('FillOval',win,[0,0,0],s.rect{1}); % small center circle, around center dot
                    Screen('FillOval',win,[255,0,0],s.rect{4}); % center dot
                end
                flipTime = Screen('Flip',win);
                flipct = flipct+1;
                data.flipTime(flipct) = flipTime;
                data.flipDate(flipct) = now;
                Screen('Close',movietexture);
            end
        else
            while (toc < blockDur/2)
                [keyIsDown, secs, keyCode, ~] = KbCheck(-3);
                if keyIsDown % If *any* key is down
                    % If t is one of the keys being pressed
                    if (ismember(KbName('t'),find(keyCode)))
                        if (secs-data.TRTime(end)) > 0.5
                            TRct = TRct+1;
                            data.TRTime(TRct) = secs;
                            data.TRDate(TRct) = now;
                            fprintf('*** T received ***\n');
                        end
                    elseif (ismember(esc,find(keyCode)))
                        Screen('PlayMovie',movie,0); % Stop Movie
                        Screen('CloseMovie',movie); % Close Move
                        clear Screen; % Clear the screen
                        Screen('CloseAll');sca; % Close the screen
                        ShowCursor;
                        ListenChar(1);
                        return
                    end
                end
                [movietexture,pts] = Screen('GetMovieImage', win, movie, 1);
                if (movietexture>0 && pts>lastpts)
                    % Store its texture handle and exact movie timestamp in
                    % arrays for later use:
                    count=count + 1;
                    data.texids(count)=movietexture;
                    data.texpts(count)=pts;
                    lastpts=pts;
                else
                    ShowCursor;
                    ListenChar(1);
                    return;
                end;
                Screen('DrawTexture',win,movietexture,[],[0 0 res.width res.height]);
                flipTime = Screen('Flip',win);
                flipct = flipct+1;
                data.flipTime(flipct) = flipTime;
                data.flipDate(flipct) = now;
                Screen('Close',movietexture);
            end
            background = [128,128,128]; % Set background color
            Screen('FillRect',win, background);
            Screen('Flip',win);
            data.blockTimes = [data.blockTimes now];
            while (toc < blockDur)
                [keyIsDown, secs, keyCode, ~] = KbCheck(-3);
                if keyIsDown % If *any* key is down
                    % If t is one of the keys being pressed
                    if (ismember(KbName('t'),find(keyCode)))
                        if (secs-data.TRTime(end)) > 0.5
                            TRct = TRct+1;
                            data.TRTime(TRct) = secs;
                            data.TRDate(TRct) = now;
                            fprintf('*** T received ***\n');
                        end
                    elseif (ismember(esc,find(keyCode)))
                        Screen('PlayMovie',movie,0); % Stop Movie
                        Screen('CloseMovie',movie); % Close Move
                        clear Screen; % Clear the screen
                        Screen('CloseAll');sca; % Close the screen
                        ShowCursor;
                        ListenChar(1);
                        return
                    end
                end
            end
        end
    end
catch ME
    Screen('PlayMovie',movie,0); % Stop Movie
    Screen('CloseMovie',movie); % Close Move
    clear Screen; % Clear the screen
    Screen('CloseAll');sca; % Close the screen
    ShowCursor;
    ListenChar;
    rethrow(ME);
end
Screen('PlayMovie',movie,0); % Stop Movie
Screen('CloseMovie',movie); % Close Move
clear Screen; % Clear the screen
Screen('CloseAll');sca; % Close the screen
ShowCursor;
ListenChar;
if ~exist(fullfile('~/Desktop/MRI',[subj '_' date]),'dir')
    mkdir(fullfile('~/Desktop/MRI',[subj '_' date]));
end
save(fullfile('~/Desktop/MRI',[subj '_' date],['run' num2str(run)]),'data','-v7.3');