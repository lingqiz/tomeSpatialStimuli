function [chex]=ringmaker_AB(blurr,wedge,res,cop)
if ~exist('blurr','var')
    blurr = 0;
end
if ~exist('wedge','var')
    wedge=30;
end
if ~exist('res','var')
    res = [1024 768];
end
% Cut off proportion
if ~exist('cop','var')
    cop = 0.01;6
end
img = zeros(2*res);
[m,n]=size(img);
[x,y]=meshgrid(1:m,1:n);
%produces exponential scales
xlin=linspace(1,n/100,n/(2*wedge));
r=fliplr(exp(xlin));
for k=1:length(r)
    distances = sqrt((x-m/2).^2+(y-n/2).^2)'<=r(k);
    angles=(radtodeg(atan2(x-m/2, y-n/2))+180)';
    for i=1:m
        for j=1:n
            %if within the distance, continue caluculation, other leave as
            %background
            if distances(i,j)==1
                if mod(floor(angles(i,j)/wedge),2)
                    img(i,j)=1; %white
                else
                    img(i,j)=-1;
                end
            end
        end
    end
    img=img*-1;
end
img(img==-1) = 0;
if blurr
    img = Blurr_AB(img,cop);
    chex = img((m/4+1):3*m/4,(n/4+1):3*n/4);
else
    chex = img((m/4+1):3*m/4,(n/4+1):3*n/4);
end
%figure;imagesc(chex);colormap gray;