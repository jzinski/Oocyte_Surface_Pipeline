function [segstack2,CH]=OocyteIDstp1V2(segstack,imstksm,tShoot)
%This identifies the largest 3-d object after binerizing and sets all
%pixels outside of it to 0.  It also sets all original 0s to 1s.


%% binarize using the smoothed images and thresh

level=double(median(imstksm(:))+round(std(double(imstksm(:)))))/256;
BW = imbinarize(imstksm,level);

% dilate the BW to fill holes
%SE = strel('cuboid',[24,24, 3]);
%BW2=imdilate(BW,SE);
%BW3=imerode(BW2,SE);

SE = strel('disk',10);
for i=1:numel(BW(1,1,:))
    curIm1=imdilate(BW(:,:,i),SE);
%    curIm1=BW(:,:,i);
    BW2(:,:,i) = imfill(curIm1,'holes');
end
%% display dilated segmeted objects
%L = bwlabeln(BW2);
%figure
%sliceViewer(L, 'colorMap', jet)

%% eliminate all except the largest object
stats = regionprops('table',BW2,'Area','PixelIdxList');
[largeObArea,largeObIdx]=max(stats.Area);

%% eliminate all except the largest object
for i=1:numel(stats.Area)
    if i~=largeObIdx
        BW2(stats.PixelIdxList{i})=0;
    end
end



%% edges
%{

for i=1:numel(imstksm(1,1,:))
    [Gmag(:,:,i),~] = imgradient(imstksm(:,:,i));
end
GmagBW=imbinarize(Gmag,0);
for i=1:numel(GmagBW(1,1,:))
  %  curIm1=imdilate(GmagBW(:,:,i),SE);
    curIm1=GmagBW(:,:,i);
    GmagBWf(:,:,i) = imfill(curIm1,'holes');
end

stats = regionprops('table',GmagBWfe,'Area','SubarrayIdx');
[largeObArea,largeObIdx]=max(stats.Area);

for i=1:numel(stats.Area)
    if i~=largeObIdx
        GmagBWf(stats.SubarrayIdx{i,1},stats.SubarrayIdx{i,2},stats.SubarrayIdx{i,3})=0;
    end
end
%}
%% watershed
%%%%%%%%%%%% smooth and binary 
h = fspecial3('ellipsoid',[30,30,round(30/8)]); %makes 3-d kernel
extrasmstk=imfilter(segstack,h,'symmetric' ); %smooth original
BWess=imbinarize(extrasmstk,round(double(median(extrasmstk(:)))/256)); %binarizae with generous boundaries
for i=1:numel(BWess(1,1,:)) %fill holes z by z
  %  curIm1=imdilate(GmagBW(:,:,i),SE);
    curIm1=BWess(:,:,i);
    BWess(:,:,i) = imfill(curIm1,'holes');
end
%%%%%%%%%%% watershed
for i=1:numel(BWess(1,1,:))
    Dess(:,:,i) = -bwdist(~BWess(:,:,i));
    Dess2(:,:,i)=imhmin(Dess(:,:,i),abs(min(min(min(Dess(:,:,i))))/7));
    Less(:,:,i)=watershed(Dess2(:,:,i),8);
end
BWcmb=BWess | BW2;
Less(BWcmb==0)=0;

%%%%%%%%% impose watershed lines onto original BW and get rid of small area
%%%%%%%%% objects
BW3=BW2 & logical(Less);

for i=1:size(BW3,3)
    BW3cur=BW3(:,:,i);
    stats = regionprops('table',BW3cur,'Centroid','Area','PixelIdxList');
    if isempty(stats)==0
    curCentroid=stats.Centroid;
    curCentroid(:,1)=curCentroid(:,1)-size(BW3,2)/2;
    curCentroid(:,2)=curCentroid(:,2)-size(BW3,1)/2;
    curDist=(curCentroid(:,1).^2+curCentroid(:,1).^2).^.5;
    areaCut=mean(stats.Area);
    curDist(stats.Area<areaCut-1)=inf;
    [~,minDIdx]=min(curDist);

    for o=1:size(stats.Centroid,1)
        if o~=minDIdx
            BW3cur(stats.PixelIdxList{o})=0;
        end
    end   
    end
    BW3(:,:,i)=BW3cur;
end

BW3=imerode(BW3,SE);

%{
BW3disp=uint8(zeros(size(BW3,1),size(BW3,2),3,size(BW3,3)));
for i=1:size(BW3,3)
    BW3disp(:,:,:,i)=imfuse(imstksm(:,:,i),BW3(:,:,i));
end
BW3disp=permute(BW3disp,[1,2,4,3]);
sliceViewer(BW3disp)
%}


%% convex hull

%convex hull
CH2=zeros(size(BW,1),size(BW,2),3,size(BW,3));
for i=1:numel(BW(1,1,:))
    CH(:,:,i) = bwconvhull(BW3(:,:,i));
    CH2(:,:,:,i)=imfuse(CH(:,:,i),segstack(:,:,i));
end
CH2=permute(CH2,[1,2,4,3]);

if tShoot==1
    figure
    sliceViewer(CH2)
end

%% set original 0s to 1s and set mask to 0s
segstack2=segstack;
segstack2(segstack2==0)=1;
segstack2(CH==0)=0;

imstksm2=imstksm;
imstksm2(imstksm2==0)=1;
imstksm2(CH==0)=0;
%{
%% resegment segstack2 using only pixels inside to determine levels
level = graythresh(imstksm2(CH~=0))*.9;
BWpostCH1= imbinarize(imstksm2,level);

level = graythresh(imstksm2(CH~=0))/2;
BWpostCH2= imbinarize(imstksm2,level);
CBWpostCHcom=[];
for i=1:numel(BW(1,1,:))
    CBWpostCHcom(:,:,:,i)=imfuse(BWpostCH1(:,:,i),BWpostCH2(:,:,i));
end
CBWpostCHcom=permute(CBWpostCHcom,[1,2,4,3]);
figure
sliceViewer(CBWpostCHcom)

%% hollow out hull

%makes pixels outside the mask NaN
imstksm2=imstksm;
imstksm2(CH==0)=NaN;


%invert and binzerize
imstksm2I=imcomplement(imstksm2);

%binerize taking new thresh each slice
for i=1:numel(imstksm2I(1,1,:))
    curIm=imstksm2I(:,:,i);
    level = graythresh(curIm(CH(:,:,i)~=0));
    level=level+(1-level)/2;
    CHBW(:,:,i) = imbinarize(curIm,level);
    levelall(i)=level;
end

%set points outside mask to 0
CHBW(CH==0)=0;

figure
sliceViewer(CHBW)

%erode then eliminate objects not attatched to center then convex hull then
%dialate
SE = strel('disk',6);
CHBW2=imerode(CHBW,SE);
sliceViewer(CHBW2)

CHBW3=imdilate(CHBW2,SE);
%}


end