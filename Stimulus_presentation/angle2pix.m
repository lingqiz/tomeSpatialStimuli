function pix = angle2pix(display,ang,dim)
%   Converts Visual Angle in degree to Pixles
%
%	Usage:
%   pix = angle2pix(display,ang,dim)
%   
%   Required:
%   display.dist (distance from screen (cm))
%   display.resolution (number of pixels of display in horizontal direction)
%   display.width (width of screen (cm)) OR display.height (height of screen (cm))
%
%   Inputs:
%   display - structure (above)
%   ang - visual angle
%
%   Optional:
%   dim - dimention, either 'height' or 'width' <default>
%
%Warning: assumes isotropic (square) pixels

%Written 11/1/07 gmb zre
% Updated 11/21/14 ASB - in case the assumption of square pixels is  false,
% and since the height is a limiting factor (i.e. max visual angle) in MRI
% stimuli, I added the 'dim' input for height or width

%% defaults
if ~exist('dim','var')
    dim = 'width';
end
%Calculate pixel size
if strcmp(dim,'width')
    pixSize = display.width/display.resolution(1);   %cm/pix
else
    pixSize = display.height/display.resolution(2);   %cm/pix
end
sz = 2*display.dist*tan(pi*ang/(2*180));  %cm
pix = round(sz/pixSize);   %pix
return

%test code
display.dist = 60; %cm
display.height = 44.5; %cm
display.resolution = [1680,1050];
ang = 2.529;
angle2pix(display,ang)




