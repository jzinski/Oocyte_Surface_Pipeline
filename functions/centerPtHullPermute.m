function [minRadBW_CH] = centerPtHullPermute(segstackM,CH,imstksm,zRatio,tShoot)
%CENTERPTHULL This plots points radiating from center
% zRatio=7.82;


%% permute inputs
SE=strel('disk',10);
CH=imdilate(CH,SE);
imstksm=permute(imstksm,[1,3,2]);
CH=permute(CH,[1,3,2]);
segstackM=permute(segstackM,[1,3,2]);

%% find center of mask at every slice
%{
for i = 1:size(CH,3)
    centroidsCur = regionprops(CH(:,:,i),'Centroid');
    if isempty(centroidsCur)==0
        centroidsAll(i,:)=[centroidsCur.Centroid];
    else
        centroidsAll(i,:)=[-1,-1];
    end
end
%}
centroidsAll = regionprops('table',CH,'Centroid'); %single centroid


%% convert each slice location to spherical coordinates using centroids

%{
%makes a meshgrid around the centroids on each slice
CH3=CH;
for i = 1:size(CH,3)
    cenCur=centroidsAll(i,:);
   [Xgrid,Ygrid]=meshgrid(1:size(segstackM, 1),1:size(segstackM, 2));
   if cenCur(1) > 0
  %     CH3(round(cenCur(2)),round(cenCur(1)),i)=0;  
 %      CH3(:,:,i)=imerode(CH3(:,:,i),SE);
        XallCart(:,:,i)=Xgrid-round(cenCur(1));
        YallCart(:,:,i)=Ygrid-round(cenCur(2));
   else
        XallCart(:,:,i)=Xgrid;
        YallCart(:,:,i)=Ygrid;
   end
end

%convert to cartesian
rhoAll=(XallCart.^2+YallCart.^2).^.5;
thetaAll=rad2deg(atan(YallCart./XallCart));
thetaAll(YallCart<=0 & XallCart>=0)=thetaAll(YallCart<=0 & XallCart>=0)+180; %this part works out quadrants to make it go from 0 to 360
thetaAll(YallCart>0 & XallCart>=0)=thetaAll(YallCart>0 & XallCart>=0)+180;
thetaAll(YallCart>0 & XallCart<0)=thetaAll(YallCart>0 & XallCart<0)+360;
%}
%makes a meshgrid around the centroids on each slice
CH3=CH;
cenCur=centroidsAll;
[Xgrid,Ygrid,Zgrid]=meshgrid(1:size(segstackM, 2),1:size(segstackM, 1),1:size(segstackM, 3));
XallCart=Xgrid-(size(CH3,2)-round(cenCur.Centroid(1)));
YallCart=Ygrid-(size(CH3,1)-round(cenCur.Centroid(2)));

% normalize z dimension
XallCart=XallCart.*zRatio;
%convert to cartesian
rhoAll=(XallCart.^2+YallCart.^2).^.5;
thetaAll=rad2deg(atan(YallCart./XallCart));
thetaAll(YallCart<=0 & XallCart>=0)=thetaAll(YallCart<=0 & XallCart>=0)+180; %this part works out quadrants to make it go from 0 to 360
thetaAll(YallCart>0 & XallCart>=0)=thetaAll(YallCart>0 & XallCart>=0)+180;
thetaAll(YallCart>0 & XallCart<0)=thetaAll(YallCart>0 & XallCart<0)+360;

%% plot all rho's that are in binary in 2D and fit spline to them
%level = graythresh(imstksm(CH==1))*1.3;
SE=strel('disk',30);
CHd=imdilate(CH,SE);
level=median(double(imstksm(CHd==1)))/256+std(double(imstksm(CHd==1)))/256;
BW = imbinarize(imstksm,level);
BW2 = bwperim(BW,8);
BW3=BW2 & CH;

curim=[];
onesRho={};
onesTheta={};
for  i=1:numel(imstksm(1,1,:))
    fitobj{i}={};
end
for i=1:numel(imstksm(1,1,:))
        curRho=rhoAll(:,:,i);
        onesRho{i}=curRho(BW3(:,:,i)==1);
        curTheta=thetaAll(:,:,i);
        onesTheta{i}=curTheta(BW3(:,:,i)==1);
        count=0;
        thetaMincur=[];
        if isempty(onesRho{i})==0
        for o=-20:380
            if o<0 || o>360
                count=count+1;
                thetaMincur(count,1)=o;
                thetaMincur(count,2)=min(onesRho{i}(round(onesTheta{i})==round(max(onesTheta{i}))));   
            else
                if isempty(min(onesRho{i}(round(onesTheta{i})==o)))==0
                    count=count+1;
                    thetaMincur(count,1)=o;
                    thetaMincur(count,2)=min(onesRho{i}(round(onesTheta{i})==o));          
                end
            end
        end
         fitobj{i}=fit(thetaMincur(:,1),thetaMincur(:,2),'smoothingspline','SmoothingParam',.001);
        end
        
end

%{
for i=1:100:numel(imstksm(1,1,:))
    if isempty(onesRho{i})==0
        figure(i)
        tiledlayout(1,2)
        nexttile
        scatter(onesTheta{i},onesRho{i})
        hold on
        axis([0,360,0,700])
        plot(fitobj{i})
        nexttile
        imshow(BW(:,:,i))
    end
end

%}
%convert the points in fitobject into a binary
minRadBW=zeros(size(segstackM, 1),size(segstackM, 2),size(segstackM, 3));
%find start and stop
objP1=min(find(max(max(CH))==1));
objP2=max(find(max(max(CH))==1));
objRange=objP2-objP1;
size1=size(CH3,1);
size2=size(CH3,2);
for z=round(objP1+objRange*1/6):round(objP2-objRange*1/6) %step z leaving the first 1/6 and last 1/6 off
    if isempty(fitobj{z})==0
    for theta=1:360 %step theta
            curRad=feval(fitobj{z},theta);  %evaluate fitted minimum
            curX=round(curRad*cosd(theta))+round(cenCur.Centroid(1)*zRatio);
            curY=round(curRad*sind(theta))+round(cenCur.Centroid(2));
            if curX<0
                curX=0;
            end
            if size2-round(curX/zRatio)>0    
               minRadBW(size1-curY,size2-round(curX/zRatio),z)=1;  %makes a plane by plane dot outline of the splines
            end
    end 
    end
    %if sind(z*10)==1
    %figure(z)
   %imshowpair(minRadBW(:,:,z),imstksm(:,:,z))
  %  end
end

%convex hull minRadBW
minRadBW_CH=zeros(size(segstackM, 1),size(segstackM, 2),size(segstackM, 3));
for z=round(objP1+objRange*1/6):round(objP2-objRange*1/6) %step z leaving the first 1/6 and last 1/6 off
    minRadBW_CH(:,:,z)=bwconvhull(minRadBW(:,:,z));  %makes a convex hull mask for the middle 2/3 of the object
end

%only does the middle 2/3 of the object


%% display output
CH2=uint8(zeros(size(BW,1),size(BW,2),3,size(BW,3)));
SE=strel('disk',2);
for i=1:numel(minRadBW_CH(1,1,:))
    fuseOut(:,:,i) = imdilate(minRadBW(:,:,i),SE);
    %CH2(:,:,:,i)=imfuse(uint8(fuseOut(:,:,i)).*2^8,imstksm(:,:,i));
    CH2(:,:,:,i)=imfuse(segstackM(:,:,i),fuseOut(:,:,i));
end
CH2=permute(CH2,[1,2,4,3]);
if tShoot==1
figure
sliceViewer(CH2)
end

%% swap dimensions back for output
minRadBW_CH=permute(minRadBW_CH,[1,3,2]);


end