function [metaOut, elePercC2,elePercC3,H1counts,H1countsC2,Cencounts,BWc2sm,BWc3sm ]=IDobjectAnglesV3(CH,c2sm,c3sm,xsize,ysize,zsize,nameStr,saveGraphs,level2,level3)
%c3 is scattered throughout the vegetal
%c2 marks animal
%old outputs:elevationC3pix,elevationCenHull2,cenC3smPixAllC3,cenC2smPixAllC3,perimHullPixAllC3,cenC2smC3,
%out
%BWc2sm=BW channel 2 bianry
%BWc3sm=BW channel 3 binary
%H1countsC2=the number of voxels in current elevation that are true for C2
%H1counts=the number of voxels in current elevation that are true for C3
%Cencounts= the total number of voxels at current elevation

%% dilate to include folicle cells

SEmine=logical(fspecial3('ellipsoid',[16,16,2]));
CHd=imdilate(logical(CH),SEmine);
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

%level(1) = graythresh(c2sm(CHd==1)); %use unfiltered stack to find cutoff
%level(2) = graythresh(c2sm);
%[~,mIdx]=max(level);

%level2=(std(double(c2sm(CHd==1)))*4+ median(double(c2sm(CHd==1))))/256;
%level2=(20+ median(double(c2sm(CHd==1))))/256;
BWc2sm = imbinarize(c2sm,level2/256);
%levelOut(1)=level2;
%level3 = graythresh(stack3(CHd==1));  %use unfiltered stack to find cutoff
%level3=(std(double(c3sm(CHd==1)))*3+ median(double(c3sm(CHd==1))))/256;
%level3=(20+ median(double(c3sm(CHd==1))))/256;
%levelOut(2)=level3;
BWc3sm = imbinarize(c3sm,level3/256);

%{
[c3HC,c3HB]=histcounts(c3sm(CHd==1),max(c3sm(CHd==1)));
tstob=fit(c3HB(1:numel(c3HC))',c3HC','poly2');
scatter(c3HB(1:numel(c3HC)),c3HC)
hold on
plot(tstob)
%}



%% eliminate signal in b2 and b3 outside of CH
BWc2smP=BWc2sm & CHd;
BWc3smP=BWc3sm & CHd;

%% get xyz of perimiter
%SE2=strel('disk',10);
SEmine2=logical(fspecial3('ellipsoid',[8,8,1]));
CHd2=imdilate(logical(CH),SEmine2);
CH2 = bwperim(CHd2,26);

%% find centroids for all
cenHull = regionprops('table',CH,'Centroid','Area','PixelList'); %single centroid
perimHull = regionprops('table',CH2,'Centroid','Area','PixelList'); %perimeter
cenC2sm = regionprops('table',BWc2smP,c2sm,'Centroid','Area','PixelList', 'PixelValues'); %single centroid
cenC3sm = regionprops('table',BWc3smP,c3sm,'Centroid','Area','PixelList','PixelValues'); %single centroid

%keep only largest c2 object
if size(cenC2sm,1)>1
    [~,maxC2]=max(cenC2sm.Area);
    cenC2sm=cenC2sm(maxC2,:);
end

%put pixels in a list
cenC3smPixAll=[];
cenC3smIntAll=[];
for i=1:numel(cenC3sm.PixelList)
    curPix=[cenC3sm.PixelList{i}];
    cenC3smPixAll=cat(1,cenC3smPixAll,curPix);
    curInts=[cenC3sm.PixelValues{i}];
    cenC3smIntAll=cat(1,cenC3smIntAll,curInts);
end

cenC2smPixAll=[];
cenC2smIntAll=[];
for i=1:numel(cenC2sm.PixelList)
    curPix=[cenC2sm.PixelList{i}];
    cenC2smPixAll=cat(1,cenC2smPixAll,curPix);
    curInts=[cenC2sm.PixelValues{i}];
    cenC2smIntAll=cat(1,cenC2smIntAll,curInts);
end

%centroid perim
perimHullPixAll=[];
for i=1:numel(perimHull.PixelList)
    curPix=[perimHull.PixelList{i}];
    perimHullPixAll=cat(1,perimHullPixAll,curPix);
end

% all C1 object volume
cenHullPixAll=[];
for i=1:numel(cenHull.PixelList)
    curPix=[cenHull.PixelList{i}];
    cenHullPixAll=cat(1,cenHullPixAll,curPix);
end

%% Make location grid with origin on centroid
[Xgrid,Ygrid,Zgrid]=meshgrid(1:size(c2sm, 2),1:size(c2sm, 1),1:size(c2sm, 3));

% correct all for center of object CH centerpoint hull
XallCart=Xgrid-size(CH,2)+round(cenHull.Centroid(1));
YallCart=Ygrid-size(CH,1)+round(cenHull.Centroid(2));
ZallCart=Zgrid-size(CH,3)+round(cenHull.Centroid(3));
cenC2smC=[cenC2sm.Centroid]-[cenHull.Centroid];
cenC3smC=[cenC3sm.Centroid]-[cenHull.Centroid];

cenC3smPixAllC=[cenC3smPixAll]-[cenHull.Centroid];
cenC2smPixAllC=[cenC2smPixAll]-[cenHull.Centroid];
perimHullPixAllC=[perimHullPixAll]-[cenHull.Centroid];
cenHullPixAllC=[cenHullPixAll]-[cenHull.Centroid];

% normalize z dimension
ZallCart=ZallCart.*zRatio;
cenC2smC(3)=cenC2smC(3).*zRatio;
cenC3smC(:,3)=cenC3smC(:,3).*zRatio;
cenC2smPixAllC(:,3)=cenC2smPixAllC(:,3).*zRatio;
cenC3smPixAllC(:,3)=cenC3smPixAllC(:,3).*zRatio;
perimHullPixAllC(:,3)=perimHullPixAllC(:,3).*zRatio;
cenHullPixAllC(:,3)=cenHullPixAllC(:,3).*zRatio;

%% Rotate so the center of the c2 centroid is north Pole
% calculate the angles of the two rotations using the center of object 2

%cenC2smC4=cenC2smC3*[1,0,0;0,cos(beta),-sin(beta);0,sin(beta),cos(beta)];
%cenC2smC4=cenC2smC3*[cos(gamma),0,sin(gamma);0,1,0;-sin(gamma),0,cos(gamma)];


%https://en.wikipedia.org/wiki/Rotation_matrix

%% Apply to centroids

theta=atan(-cenC2smC(1)/cenC2smC(2));
cenC2smC2=cenC2smC*[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1];
beta=atan(cenC2smC(3)/cenC2smC2(2))+pi/2;
%beta=atan(cenC2smC(3)/-(cenC2smC(1)^2+cenC2smC(2)^2)^.5)+pi/2;
cenC2smC3=cenC2smC2*[1,0,0;0,cos(beta),-sin(beta);0,sin(beta),cos(beta)];
flipYes=0;
if cenC2smC3(3)<0
    flipYes=1;
end
if flipYes==1
        cenC2smC3(:,3)=cenC2smC3(:,3)*-1;
end

rot1=[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1];
rot2=[1,0,0;0,cos(beta),-sin(beta);0,sin(beta),cos(beta)];

cenC3smC2=cenC3smC*rot1;
cenC3smC3=cenC3smC2*rot2;
if flipYes==1
    cenC3smC3(:,3)=cenC2smC3(:,3)*-1;
end

%{
figure
scatter3(cenC2smC3(:,1),cenC2smC3(:,2),cenC2smC3(:,3))
hold on
scatter3(cenC3smC3(:,1),cenC3smC3(:,2),cenC3smC3(:,3),80,[1,0,0],'.')
axis([min(min(min(XallCartProt))),max(max(max(XallCartProt))),min(min(min(YallCartProt))),max(max(max(YallCartProt))),min(min(min(ZallCartProt))),max(max(max(ZallCartProt)))])
axis equal
%}
%{
%% apply to grids
allCart=[XallCart(:),YallCart(:),ZallCart(:)];

allCartRot1=allCart*rot1;
allCartRot1=allCartRot1*rot2;

XallCartProt = reshape(allCartRot1(:,1), size(XallCart,1), size(XallCart,2),size(XallCart,3));
YallCartProt = reshape(allCartRot1(:,2), size(XallCart,1), size(XallCart,2),size(XallCart,3));
ZallCartProt = reshape(allCartRot1(:,3), size(XallCart,1), size(XallCart,2),size(XallCart,3));



%% convert to spherical
rhoAll=(XallCartProt.^2+YallCartProt.^2+ZallCartProt.^2).^.5;
thetaAll=rad2deg(atan(YallCartProt./XallCartProt));
thetaAll(YallCartProt<=0 & XallCartProt>=0)=thetaAll(YallCartProt<=0 & XallCartProt>=0)+180; %this part works out quadrants to make it go from 0 to 360
thetaAll(YallCartProt>0 & XallCartProt>=0)=thetaAll(YallCartProt>0 & XallCartProt>=0)+180;
thetaAll(YallCartProt>0 & XallCartProt<0)=thetaAll(YallCartProt>0 & XallCartProt<0)+360;
phiAll=rad2deg(acos(ZallCartProt./(rhoAll)));

%% Plot
histogram(phiAll(BWc2smP==1))
histogram(rhoAll(BWc2smP==1))
histogram(thetaAll(BWc2smP==1))
%}
%% Apply rotation to pixels
%c2 pixels
cenC3smPixAllC2=cenC3smPixAllC*[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1]; %rotate channel 2 object pixels
cenC3smPixAllC3=cenC3smPixAllC2*[1,0,0;0,cos(beta),-sin(beta);0,sin(beta),cos(beta)];
if flipYes==1
    cenC3smPixAllC3(:,3)=cenC3smPixAllC3(:,3)*-1;
end

% c3 pixels
cenC3smPixAllC2=cenC3smPixAllC*[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1]; %rotate channel 2 object pixels
cenC3smPixAllC3=cenC3smPixAllC2*[1,0,0;0,cos(beta),-sin(beta);0,sin(beta),cos(beta)];
if flipYes==1
    cenC3smPixAllC3(:,3)=cenC3smPixAllC3(:,3)*-1;
end

cenC2smPixAllC2=cenC2smPixAllC*[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1]; %rotate channel 3 object pixels
cenC2smPixAllC3=cenC2smPixAllC2*[1,0,0;0,cos(beta),-sin(beta);0,sin(beta),cos(beta)];
if flipYes==1
    cenC2smPixAllC3(:,3)=cenC2smPixAllC3(:,3)*-1;
end

perimHullPixAllC2=perimHullPixAllC*[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1]; %rotate perim object pixels
perimHullPixAllC3=perimHullPixAllC2*[1,0,0;0,cos(beta),-sin(beta);0,sin(beta),cos(beta)];
if flipYes==1
    perimHullPixAllC3(:,3)=perimHullPixAllC3(:,3)*-1;
end

cenHullPixAllC2=cenHullPixAllC*[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1]; %rotate all CG object pixels
cenHullPixAllC3=cenHullPixAllC2*[1,0,0;0,cos(beta),-sin(beta);0,sin(beta),cos(beta)];
if flipYes==1
    cenHullPixAllC3(:,3)=cenHullPixAllC3(:,3)*-1;
end
%% Convert to spherical
%{
rhoAllCen=(cenC3smC3(:,1).^2+cenC3smC3(:,2).^2+cenC3smC3(:,3).^2).^.5;
thetaAllCen=rad2deg(atan(cenC3smC3(:,2)./cenC3smC3(:,1)));
thetaAllCen(cenC3smC3(:,2)<=0 & cenC3smC3(:,1)>=0)=thetaAllCen(cenC3smC3(:,2)<=0 & cenC3smC3(:,1)>=0)+180; %this part works out quadrants to make it go from 0 to 360
thetaAllCen(cenC3smC3(:,2)>0 & cenC3smC3(:,1)>=0)=thetaAllCen(cenC3smC3(:,2)>0 & cenC3smC3(:,1)>=0)+180;
thetaAllCen(cenC3smC3(:,2)>0 & cenC3smC3(:,1)<0)=thetaAllCen(cenC3smC3(:,2)>0 & cenC3smC3(:,1)<0)+360;
phiAllCen=rad2deg(acos(cenC3smC3(:,3)./(rhoAllCen)));

rhoAllCen2=(cenC2smC3(:,1).^2+cenC2smC3(:,2).^2+cenC2smC3(:,3).^2).^.5;
thetaAllCen2=rad2deg(atan(cenC2smC3(:,2)./cenC2smC3(:,1)));
thetaAllCen2(cenC2smC3(:,2)<=0 & cenC2smC3(:,1)>=0)=thetaAllCen(cenC2smC3(:,2)<=0 & cenC2smC3(:,1)>=0)+180; %this part works out quadrants to make it go from 0 to 360
thetaAllCen2(cenC2smC3(:,2)>0 & cenC2smC3(:,1)>=0)=thetaAllCen(cenC2smC3(:,2)>0 & cenC2smC3(:,1)>=0)+180;
thetaAllCen2(cenC2smC3(:,2)>0 & cenC2smC3(:,1)<0)=thetaAllCen(cenC2smC3(:,2)>0 & cenC2smC3(:,1)<0)+360;
phiAllCen2=rad2deg(acos(cenC2smC3(:,3)./(rhoAllCen2)));
%}
%centroids c2
[azimuthC2r,elevationC2r,rC2r] = cart2sph(cenC2smC3(:,1),cenC2smC3(:,2),cenC2smC3(:,3));
elevationC2=-rad2deg(elevationC2r)+90;
azimuthC2=rad2deg(azimuthC2r);

%pixels hull
[azimuthCenHull,elevationCenHull,rCenHull] = cart2sph(cenHullPixAllC3(:,1),cenHullPixAllC3(:,2),cenHullPixAllC3(:,3));
elevationCenHull2=-rad2deg(elevationCenHull)+90;
azimuthCenHull2=rad2deg(azimuthCenHull);

%perim hull
[azimuthPerimHull,elevationPerimHull,rPerimHull] = cart2sph(perimHullPixAllC3(:,1),perimHullPixAllC3(:,2),perimHullPixAllC3(:,3));
elevationPerimHull2=-rad2deg(elevationPerimHull)+90;
azimuthPerimHull2=rad2deg(azimuthPerimHull);

%{
[azimuthC3r,elevationC3r,rC3r] = cart2sph(cenC3smC3(:,1),cenC3smC3(:,2),cenC3smC3(:,3));
elevationC3=rad2deg(elevationC3r)-90;
azimuthC3=rad2deg(azimuthC3r);
%}
%pixelc2
[azimuthC2pix,elevationC2pix,rC2pix] = cart2sph(cenC2smPixAllC3(:,1),cenC2smPixAllC3(:,2),cenC2smPixAllC3(:,3));
elevationC2pix=-rad2deg(elevationC2pix)+90;
azimuthC2pix=rad2deg(azimuthC2pix);

%pixelc3
[azimuthC3pix,elevationC3pix,rC3pix] = cart2sph(cenC3smPixAllC3(:,1),cenC3smPixAllC3(:,2),cenC3smPixAllC3(:,3));
elevationC3pix=-rad2deg(elevationC3pix)+90;
azimuthC3pix=rad2deg(azimuthC3pix);


%% add meta data to xls
% add volume of C1, volume of C2, Volume of C3
%means are summaries of pixel values
%volumes are converted from x=voxels using xsize, ysize, zsize (um^3)
%metaOut=[volumeC1,volumeC2,volumeC3,meanIntC2,meanIntC3];

volumeC1=sum(sum(sum(CHd)))*xsize*ysize*zsize;
volumeC2=sum(sum(sum(size(cenC2smPixAll,1))))*xsize*ysize*zsize;
volumeC3=sum(sum(sum(size(cenC3smPixAll,1))))*xsize*ysize*zsize;
meanIntC2=mean(cenC2smIntAll);
meanIntC3=mean(cenC3smIntAll);
metaOut=[volumeC1,volumeC2,volumeC3,meanIntC2,meanIntC3,level2,level3];



%{
% centroids pre-rotate
figure
scatter3(cenC2smC(:,1),cenC2smC(:,2),cenC2smC(:,3),20,[0,0,1],'.')
hold on
scatter3(cenC3smC(:,1),cenC3smC(:,2),cenC3smC(:,3),20,[1,0,0],'.')
%plot3([0,],[0,],[0,],'-')
axis equal

% centroids post-rotate
figure
scatter3(cenC2smC3(:,1),cenC2smC3(:,2),cenC2smC3(:,3),20,[0,0,1],'.')
hold on
scatter3(cenC3smC3(:,1),cenC3smC3(:,2),cenC3smC3(:,3),20,[1,0,0],'.')
%plot3([0,],[0,],[0,],'-')
axis equal
%}

%{
figure
histogram(azimuthC3pix,180)
grid on;
xlim([-180, 180]);
xlabel('Azimuth (degrees)', 'FontSize', 14);
ylabel('Pixel Count', 'FontSize', 14);
%}

%% gif of c2

if saveGraphs==1
    c2smS=uint8(double(c2sm).*256/double(max(max(max(c2sm)))));
    CH4=uint8(zeros(size(BWc2smP,1),size(BWc2smP,2),3,size(BWc2smP,3)));
    for i=1:numel(BWc2smP(1,1,:))
        CH4(:,:,:,i)=imfuse(c2smS(:,:,i),uint8(BWc2smP(:,:,i)).*70, 'Scaling', 'none');
    end
    CH4=permute(CH4,[1,2,4,3]);
    figure
    s=sliceViewer(CH4);
    
    hAx = getAxesHandle(s);

    filename = [nameStr '_c2.gif'];

    sliceNums = 1:numel(BWc2smP(1,1,:));

    for idx = sliceNums
        % Update slice number
        s.SliceNumber = idx;
        % Use getframe to capture image
        I = getframe(hAx);
        [indI,cm] = rgb2ind(I.cdata,256);
        % Write frame to the GIF file
        if idx == 1
            imwrite(indI,cm,filename,'gif','Loopcount',inf,'DelayTime', 0.1);
        else
            imwrite(indI,cm,filename,'gif','WriteMode','append','DelayTime', 0.1);
        end
    end
    
end

%% gif of c3

if saveGraphs==1
    c3smS=uint8(double(c3sm).*256/double(max(max(max(c3sm)))));
    CH4=uint8(zeros(size(BWc2smP,1),size(BWc2smP,2),3,size(BWc2smP,3)));
    for i=1:numel(BWc2smP(1,1,:))
        CH4(:,:,:,i)=imfuse(c3smS(:,:,i),uint8(BWc3smP(:,:,i)).*70, 'Scaling', 'none');
    end
    CH4=permute(CH4,[1,2,4,3]);
    figure
    s=sliceViewer(CH4);
    
    hAx = getAxesHandle(s);

    filename = [nameStr '_c3.gif'];

    sliceNums = 1:numel(BWc2smP(1,1,:));

    for idx = sliceNums
        % Update slice number
        s.SliceNumber = idx;
        % Use getframe to capture image
        I = getframe(hAx);
        [indI,cm] = rgb2ind(I.cdata,256);
        % Write frame to the GIF file
        if idx == 1
            imwrite(indI,cm,filename,'gif','Loopcount',inf,'DelayTime', 0.1);
        else
            imwrite(indI,cm,filename,'gif','WriteMode','append','DelayTime', 0.1);
        end
    end
    
end



%% determine how much of surface is covered as a binary

%plot c3 
C3perimCoverage=[];
count=0;
for i=0:5:175
    for o=-180:5:175
        count=count+1;
        curNumR=find(elevationC3pix>=i & elevationC3pix<i+5 & azimuthC3pix>=o & azimuthC3pix<o+5);
        C3perimCoverage(count,:)=[i,o,~isempty(curNumR)];
    end
end
figure
subplot(3,3,3)
scatter(C3perimCoverage(:,1),C3perimCoverage(:,2),100,C3perimCoverage(:,3),'.')
colormap winter
ylim([-180,180])
xlim([0,180])
xlabel('Angle From Animal Pole (degrees)', 'FontSize', 12);
ylabel('Angle Around Margin (degrees)', 'FontSize', 12);
title('Dazl Coverage')

%set(gca, 'Xdir','reverse')
count=0;
for i=0:5:175
    count=count+1;
    elePercC3(count)=sum(C3perimCoverage(i<=C3perimCoverage(:,1) & C3perimCoverage(:,1)<i+5,3))/numel(C3perimCoverage(i<=C3perimCoverage(:,1) & C3perimCoverage(:,1)<i+5,3));
end
subplot(3,3,6)
plot(0:5:175,elePercC3,'LineWidth', 3, 'Color', [.75,0,0])
xlim([0,180])
xlabel('Angle From Animal Pole (degrees)', 'FontSize', 12);
ylabel('Fraction Surface Covered by Dazl', 'FontSize', 12);
title(['Dazl Fraction Coverage=' num2str(sum(elePercC3)/numel((elePercC3)),2)])


% pixels c3 histogram normalized by volume
subplot(3,3,9)
%[Perimcounts,Perimedges]=histcounts(elevationPerimHull2,[0:5:180]);
[Cencounts,Cenedges]=histcounts(elevationCenHull2,[0:5:180]);  %area of all pixels
[H1counts,H1edges]=histcounts(elevationC3pix,[0:5:180]);   %channel 3 pixels
%H1volume=(H1counts)*xsize*ysize*zsize.*(max(Cencounts)./Cencounts);  %normalizing for total amount of pixels
H1volume=(H1counts)./Cencounts;  %normalizing for total amount of pixels
%H1volume=(H1counts)*xsize*ysize*zsize.*(max(Perimcounts)./Perimcounts);
bar(H1edges(1:end-1),H1volume, 'r')
hold on
grid on;
xlim([0, 180]);
xlabel('Angle From Animal Pole (degrees)', 'FontSize', 12);
ylabel('Percent Volume C3', 'FontSize', 12);
title(['C3 Volume= ' num2str(round(volumeC3)) 'um^3'] ,  'Interpreter', 'none')

%plot c2
C2perimCoverage=[];
count=0;
for i=0:5:175
    for o=-180:5:175
        count=count+1;
        curNumR=find(elevationC2pix>=i & elevationC2pix<i+5 & azimuthC2pix>=o & azimuthC2pix<o+5);
        C2perimCoverage(count,:)=[i,o,~isempty(curNumR)];
    end
end
subplot(3,3,2)
scatter(C2perimCoverage(:,1),C2perimCoverage(:,2),100,C2perimCoverage(:,3),'.')
%set(gca, 'XDir','reverse')
ylim([-180,180])
xlim([0,180])
colormap winter
xlabel('Angle From Animal Pole (degrees)', 'FontSize', 12);
ylabel('Angle Around Margin (degrees)', 'FontSize', 12);
title('Cyclin B1 Coverage')

count=0;
for i=0:5:175
    count=count+1;
    elePercC2(count)=sum(C2perimCoverage(i<=C2perimCoverage(:,1) & C2perimCoverage(:,1)<i+5,3))/numel(C2perimCoverage(i<=C2perimCoverage(:,1) & C2perimCoverage(:,1)<i+5,3));
end
subplot(3,3,5)
plot(0:5:175,elePercC2, 'LineWidth', 3, 'Color', [0,.75,0])
xlim([0,180])
xlabel('Angle From Animal Pole (degrees)', 'FontSize', 12);
ylabel('Fraction Surface Covered by CyclinB1', 'FontSize', 12);
title(['Cyclin B1 Fraction Coverage=' num2str(sum(elePercC2)/numel((elePercC2)),2)])


% pixels c2 histogram normalized by volume
subplot(3,3,8)
%[Perimcounts,Perimedges]=histcounts(elevationPerimHull2,[0:5:180]);
[H1countsC2,H1edgesC2]=histcounts(elevationC2pix,[0:5:180]);   %channel 3 pixels
%H1volume=(H1counts)*xsize*ysize*zsize.*(max(Cencounts)./Cencounts);  %normalizing for total amount of pixels
H1volumeC2=(H1countsC2)./Cencounts;  %normalizing for total amount of pixels
%H1volume=(H1counts)*xsize*ysize*zsize.*(max(Perimcounts)./Perimcounts);
bar(H1edges(1:end-1),H1volumeC2, 'g')
hold on
grid on;
xlim([0, 180]);
xlabel('Angle From Animal Pole (degrees)', 'FontSize', 12);
ylabel('Percent Volume C2', 'FontSize', 12);
title(['C2 Volume= ' num2str(round(volumeC2)) 'um^3'] ,  'Interpreter', 'none')

%% plot
% 3d pixels pre-rotate
subplot(3,3,4)
scatter3(cenC3smPixAllC(:,1),cenC3smPixAllC(:,2),cenC3smPixAllC(:,3),20,[.75,0,0],'.')
hold on
scatter3(cenC2smPixAllC(:,1),cenC2smPixAllC(:,2),cenC2smPixAllC(:,3),20,[0,.75,0],'.')
scatter3(perimHullPixAllC(1:500:end,1),perimHullPixAllC(1:500:end,2),perimHullPixAllC(1:500:end,3),20,[.8,.8,.8],'.')
plot3([cenC2smC(:,1),0,-cenC2smC(:,1)],[cenC2smC(:,2),0,-cenC2smC(:,2)],[cenC2smC(:,3),0,-cenC2smC(:,3)],'-','LineWidth',4,'Color', [0,0,0])
axis equal
title(['C2 thresh=' num2str(round(level2)) '; C3 thresh=' num2str(round(level3))])
view(0,-90)


% 3d pixels post-rotate

subplot(3,3,7)
scatter3(cenC3smPixAllC3(:,1),cenC3smPixAllC3(:,2),cenC3smPixAllC3(:,3),20,[.75,0,0],'.')
hold on
scatter3(cenC2smPixAllC3(:,1),cenC2smPixAllC3(:,2),cenC2smPixAllC3(:,3),20,[0,.75,0],'.')
scatter3(perimHullPixAllC3(1:500:end,1),perimHullPixAllC3(1:500:end,2),perimHullPixAllC3(1:500:end,3),20,[.8,.8,.8],'.')
%scatter3(perimHullPixAllC3(1:end,1),perimHullPixAllC3(1:end,2),perimHullPixAllC3(1:end,3),20,[.8,.8,.8],'.')
plot3([cenC2smC3(:,1),0,-cenC2smC3(:,1)],[cenC2smC3(:,2),0,-cenC2smC3(:,2)],[cenC2smC3(:,3),0,-cenC2smC3(:,3)],'-','LineWidth',4,'Color', [0,0,0])
%scatter3(0,0,0,500,[1,0,0],'.')
axis equal
caxis([0 180]);
view(0,0)

title(['C1 Volume= ' num2str(round(volumeC1)) 'um^3'] ,  'Interpreter', 'none')
subplot(3,3,1)
%{
c3smS=uint8(double(c3sm).*256/double(max(max(max(c3sm)))));
maxProjC3=uint8(zeros(size(BWc2smP,1),size(BWc2smP,2),3,size(BWc2smP,3)));
for i=1:numel(BWc2smP(1,1,:))
    maxProjC3(:,:,:,i)=imfuse(c3smS(:,:,i),uint8(BWc3smP(:,:,i)).*70, 'Scaling', 'none');
end
imagesc(max(maxProjC3,[],4))
%}
maxProjAll(:,:,1)=uint8(rescale(max(c3sm,[],3), 0,256,'InputMax',40));
maxProjAll(:,:,2)=uint8(rescale(max(c2sm,[],3), 0,256,'InputMax',60));
maxProjAll(:,:,3)=uint8(zeros(size(c2sm(:,:,1))));
imagesc(maxProjAll)
title(nameStr,  'Interpreter', 'none')
set(gcf,'units','points','position',[10,10,900,900])

if saveGraphs==1
saveas(gcf, [nameStr '_masterplot_all.tif'])
end

end