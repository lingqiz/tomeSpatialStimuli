function [chex]=checker(show,display,rings,wedges,radmult,square,cpd)
%%   Creates a radial checkerboard
%
%   [chex]=checker(show,display,rings,wedges,radmult,square,cpd)
%
%   inputs:
%       show = show image (default = 0; don't show)
%       display = structure containing the following parameters
%           display.dist = distance from screen (cm)
%           display.width = width of screen (cm)
%           display.resolution = resolution of screen (pixels)
%       rings = number of rings
%       wedges = number of wedges
%       radmult = radian multiplier (needed to expand the center image) 
%       square = 0=radial <default> or 1=square checkerboard
%       cpd = cycles per degree, 2 cpd <default>
%   outputs:
%       chex = checkerboard image
%
%   Usage:
%       display.dist = 124.5;
%       display.width = 50.4;
%       display.resolution = [1024 768];
%       chex = checker(1,display,10,20,.205,0,.5)
%
%   Written by Andrew S Bock July 2014

%%
if ~exist('show','var')
    show = 0;
end
if ~exist('display','var')
    display.dist = 124.25; % distance from screen (cm)
    display.width = 50.4; % width of screen (cm)
    display.resolution = [1024 768];
end
if ~exist('rings','var')
    rings = 30;
end
if ~exist('wedges','var')
    wedges = 2*rings;
end
if ~exist('radmult','var')
    radmult = .25;
end
if ~exist('square','var')
    square = 0;
end
if ~exist('cpd','var')
    cpd = 0.5;
end
%% Make chex
X = display.resolution(1);
Y = display.resolution(2);
[xpx,ypx] = meshgrid((1:X)-X/2,(1:Y)-Y/2);
x = pix2angle(display,xpx);
y = pix2angle(display,ypx);
[ang,rad] = cart2pol(x,y);
chex = sign(sin((wedges/2)*ang).*sin(((rings+1)*3)./(exp(radmult*rad))));
if square
chex = sign(cos(2*pi*x*cpd) .* cos(2*pi*y*cpd)); 
end
chex(chex==-1) = 0;
if show
    figure;imagesc(chex);colormap gray
end