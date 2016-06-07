function [ SGmask ] = SGaussMask( SzRows, SzCols, Width, Invert)
%
% for inverting, either nothing or N is false, anything else is true;
if nargin < 4
    Invert=false;
else
    if strcmpi(Invert,'N')
        Invert = false;
    else
        Invert = true;
    end
end
SzRows=round(abs(SzRows));
SzCols=round(abs(SzCols));
if nargin < 3
    Width=0.3*min(SzRows, SzCols);
end
Width=double(Width);

SGmask=double(zeros(SzRows, SzCols));
centrow=double(SzRows)/2;
centcol=double(SzCols)/2;
for i=1:SzRows
    for j=1:SzCols
        SGmask(i,j)=(exp(-((i-centrow)^2+(j-centcol)^2)/(2*Width^2)));% / (2*pi*Width^2);
        %SGmask(i,j)=exp(-((i-centrow)^2+(j-centcol)^2)/2*Width^2) / Width*sqrt(2*pi);
%         if sqrt((i-centrow)^2+(j-centcol)^2)<Width
%             SGmask(i,j) = 1;
%         end
    end
end
if Invert
    SGmask=double(ones(SzRows, SzCols))-SGmask;
end
