function [lpimg,lpfilt]=blur(show,display,img,cpd,pad)
%%   Blurs and image using a low-pass filter
%
%   [lpimg,lpfilt] = blur(show,display,img,cpd,pad)
%
%   inputs:
%       show = show image (default = 0; don't show)
%       display = structure containing the following parameters
%           display.dist = distance from screen (cm)
%           display.width = width of screen (cm)
%           display.resolution = resolution of screen (pixels)
%       img = 2D image matrix;
%       cpd = cycles/degree cutoff
%       pad = pad image by mirroring edges, 0 = no pad <default>
%   outputs:
%       lpimg = low-pass filtered image
%   Usage:
%       lpimg = blur(1,display,img,0.01,0)
%
%   Written by Andrew S Bock July 2014
%% Cut off proportion
if ~exist('pad','var')
    pad = 0;
end
if ~exist('display','var')
    display.dist = 124.25; % distance from screen (cm)
    display.width = 50.4; % width of screen (cm)
    display.resolution = [1024 768];
end
if ~exist('cpd','var')
    cpd = 0.5;
end
if ~exist('show','var')
    show = 0;
end
%% Calculate size of image, and cut off frequency based on cycles/deg
[H,W] = size(img);
cutfreq = round((H/angle2pix(display,1))*cpd+1);
%% Pad image
if pad
    p = round(0.25*W);
    pimg = zeros(H+2*p,W+2*p);
    pimg(p+1:p+H,p+1:p+W) = img;
    %Top and Bottom
    pimg(1:p, p+1:p+W) = repmat(img(1,1:end), p, 1);
    pimg(p+H+1:end,p+1:p+W) = repmat(img(end,1:end), p, 1);
    %Left and Right
    pimg(p+1:p+H, 1:p) = repmat(img(1:end,1), 1, p);
    pimg(p+1:p+H, p+W+1:end) = repmat(img(1:end,end), 1, p);
    %Corners
    pimg(1:p, 1:p) = img(1,1); %Top-left
    pimg(1:p, p+W+1:end) = img(1,end); %Top-right
    pimg(p+H+1:end, 1:p) = img(end,1); %Bottom-left
    pimg(p+H+1:end,p+W+1:end) = img(end,end); %Bottom-right
else
    pimg = img;
end
%% Filter image
fftA=fft2(pimg);
sfftA=fftshift(fftA);
[rows, cols]=size(sfftA); %get the number of rows and columns.
%Create the masks for filters
lpfilt=SGaussMask(rows,cols,cutfreq,'N');
lpfftA=ifftshift(sfftA.*lpfilt);
lpimg=abs(ifft2(lpfftA));
%% Remove padding
if pad
    lpimg = lpimg(p+1:p+H,p+1:p+W);
    lpfilt = lpfilt(p+1:p+H,p+1:p+W);
end
if show
    figure;imagesc(lpimg);colormap gray
end



