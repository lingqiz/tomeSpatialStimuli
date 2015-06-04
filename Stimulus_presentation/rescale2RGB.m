function newimg = rescale2RGB(img,pContrast)
% rescale img (-1 to 1) to RGB (0 to 255)
newimg = round((pContrast*img+1)*127.5); 