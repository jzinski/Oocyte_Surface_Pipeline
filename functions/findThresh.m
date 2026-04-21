function [medianSm,meanSm,otsuSm,stdSm]=findThresh(CH,c2sm,c3sm,xsize,ysize,zsize)
%c3 is scattered throughout the vegetal
%c2 marks animal

%% dilate to include folicle cells
%SE=strel('disk',20);
SE=logical(fspecial3('ellipsoid',[16,16,2]));
CHd=imdilate(logical(CH),SE);
zRatio=zsize/xsize;

%% use k-means to find boundaries
%{
c3smCH=c3sm;
c3smCH(CHd==0)=median(c3sm(:));
L3 = imsegkmeans3(c3smCH,3);
BWc3sm=L3==3;

c2smCH=c2sm;
c2smCH(CHd==0)=median(c2sm(:));
L2 = imsegkmeans3(c2smCH,3);
BWc2sm=L2==3;
%}

%% find binaries of c2 and c3 using otsu


medianSm(1)=median(c2sm(CHd==1));
meanSm(1)=mean(c2sm(CHd==1));
otsuSm(1)= graythresh(c2sm(CHd==1));
stdSm(1)= std(double(c2sm(CHd==1)));

medianSm(2)=median(c3sm(CHd==1));
meanSm(2)=mean(c3sm(CHd==1));
otsuSm(2)= graythresh(c3sm(CHd==1));
stdSm(2)= std(double(c3sm(CHd==1)));


end