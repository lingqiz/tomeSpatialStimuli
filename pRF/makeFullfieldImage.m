function img = makeFullfieldImage( display, stim )
% img = makeFullfieldImage( display, stim )
%
% DESCRIPTION
% returns the basic image of the stimulus to be used
% INPUT
%   display
%       .screenAngle - resolution of screen in degrees
%       .resolution  - resolution of screen in pixels
%       .squares - 'square' or 'ecc' for standard or radial checkerboard
%   center         - bar center - default is 0
% OUTPUT
%   img - stimulus image [display.resolution] ranging between 0 and 1
%
% 04/2013 PB

if ~exist( 'center', 'var' )
    center = [0 0];
end

if ~isfield(stim,'ang')
    stim.ang = 18;
end

if ~isfield(stim,'rad')
    stim.rad = 4;
end

if ~isfield(stim,'squares')
    stim.squares = 'square'; % default to standard square checkerboard
end

% make basis functions
[xpx ypx] = meshgrid(   ( 1 : display.resolution(1) ) - display.resolution(1)/2 ,...
    ( 1 : display.resolution(2) )  - display.resolution(2)/2 );
x = pix2angle( display, xpx );
y = pix2angle( display, ypx );

[ang,rad] = cart2pol(x,y);

%mask = ang < -pi/2 | ang > pi/2;

% add checkerboard pattern
if strcmp(stim.squares,'square')
    chex = sign( cos(2*pi*x*stim.cyclePerDeg) .* cos(2*pi*y*stim.cyclePerDeg) ); %Standard checkerboard
elseif strcmp(stim.squares,'ecc')
    chex = sign( sin(stim.ang*ang) .* sin(stim.rad*rad)); %Size of 'squares' increase with eccentcity
end

img = chex;
%img = chex .* mask;
