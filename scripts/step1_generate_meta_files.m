%% Path setup - adds ../functions and ../thirdparty to MATLAB path
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir,'..','functions'));
addpath(genpath(fullfile(scriptDir,'..','thirdparty')));

%% Data folder - set to the folder containing your .tif stacks and metadata xlsx
dataDir = fullfile(scriptDir,'..','data');
cd(dataDir);

tableName='Oocyte_sample_list';
metaAll=readtable([tableName '.xlsx']);
%to do: put a convex full step in combine binaries?? 

%% gaussian filter settings (keep same as master function)
x1=24;
y1=24;
z1=3;

x2=24;
y2=24;
z2=3;

x3=14;
y3=14;
z3=3;
tShoot=0;
%% x y z pixel size microns (user can change)
xsize=0.3321;
ysize=0.3321;
zsize=2.6;
zRatio=zsize/xsize;

%% display option
% 1 is yes, 0 is no
individOocyteOutlines=1; %display the segmented oocyte stack
tShoot=0; %extra images for troubleshooting segmentation
saveGraphs=1; %save the histgrams to file

%% determine directory data
dirList=dir;
names_cell = {dirList.name};
tableFinal=[];

%% for each unique date find a thresh
allUni=unique(metaAll.date);
for i=1:numel(allUni)
    curDate=allUni(i);
    dateIdx=metaAll.date==curDate; %index of rows in this set
    
    curTable=metaAll(dateIdx,:);
    
    %prebuild the stores for output variables
    medianAll=zeros(numel(curTable.name),2);
    meanAll=zeros(numel(curTable.name),2);
    otsuAll=zeros(numel(curTable.name),2);
    stdAll=zeros(numel(curTable.name),2);
    passAll=zeros(numel(curTable.name),1);
    
    %% find histogram of c2 and c3 pixels in this group
    
    for u=1:numel(curTable.name) %loops through all names from the current date and generates a variables which will be added to the current table
        %clear old variables
        clear imstksm segstack minRadBW_CH_Tot
        skipMeans=0;
        
        %import
        for o=1:3
            [segstack(:,:,:,o)]=importstackone(curTable.name{u}, 3, o); %imports current oocyte
        end
        %smooth
        [imstksm(:,:,:,1)]=smooth3D(segstack(:,:,:,1), x1, y1, z1);  %Apply Gaussian Filter c1
        [imstksm(:,:,:,2)]=smooth3D(segstack(:,:,:,2), x2, y2, z2);  %Apply Gaussian Filter c2
        [imstksm(:,:,:,3)]=smooth3D(segstack(:,:,:,3), x3, y3, z3);  %Apply Gaussian Filter c3

        %% identify oocyte
        if any(strcmp(names_cell,[curTable.name{u} '_c1mask.mat']))  %checks if a .mat c1 binary 3d mask exists for this files
            % loads binary if it exists
            load([curTable.name{u} '_c1mask.mat'])
            disp('Oocyte ID Load Complete')
        else
            try
            % performs binary generation if it does not exist
            [segstackM,CH]=OocyteIDstp1V2(segstack(:,:,:,1),imstksm(:,:,:,1),tShoot); %identify oocyte
            [minRadBW_CH,minRadBW_CH_full] = centerPtHullV3(segstackM(:,:,:,1),CH,imstksm(:,:,:,1),tShoot);
            CH2=logical(minRadBW_CH_full) & CH;
            [minRadBW_CHp] = centerPtHullPermute(segstackM,CH2,imstksm(:,:,:,1),zRatio,tShoot); %makes a 3-D binary of missing pieces by turning it on its side and repeating
            minRadBW_CHp=minRadBW_CHp & CH2;
            [minRadBW_CH_Tot]=combineBinaries(minRadBW_CH, minRadBW_CHp,imstksm(:,:,:,1),individOocyteOutlines,curTable.name{u}); %combines the two binaries
            disp('Oocyte ID Complete')
            catch
                skipMeans=1;   %sets to skip saving the mask and taking means
                disp('Oocyte ID Abort')
            end
        end
        %% identidy the common c2 and c3 threshold based on pixel values in identified embryos
        if skipMeans==0 
            %get means and save mask
            [medianSm,meanSm,otsuSm,stdSm]=findThresh(minRadBW_CH_Tot,imstksm(:,:,:,2),imstksm(:,:,:,3),xsize,ysize,zsize);         
            save([curTable.name{u} '_c1mask.mat'],'minRadBW_CH_Tot')  % save mask
                            %add in something to extract regionprops3 but accounting
                %for the differences in scale between xy and z
             %   C1_props=regionprops3(minRadBW_CH_Tot);
            passCur=1;
        else
            %skip save mask and set all measurements to 0
            medianSm=[0,0];
            meanSm=[0,0];
            otsuSm=[0,0];
            stdSm=[0,0];
            passCur=0;
        end
        %% log values for variables
        medianAll(u,:)=medianSm; %col 1 is c2, col 2 is c3
        meanAll(u,:)=meanSm; %col 1 is c2, col 2 is c3
        otsuAll(u,:)=otsuSm; %col 1 is c2, col 2 is c3 
        stdAll(u,:)=stdSm; %col 1 is c2, col 2 is c3
        passAll(u)=passCur;

    end
    %% enter new values into table  (calibrated for 8-bit)
    name = 'mean_C2';
    curTable.(name) =meanAll(:,1);
    name = 'mean_C3';
    curTable.(name) =meanAll(:,2);
    name = 'median_C2';
    curTable.(name) =medianAll(:,1);
    name = 'median_C3';
    curTable.(name) =medianAll(:,2);
    name = 'otsu_C2';
    curTable.(name) =otsuAll(:,1).*256;
    name = 'otsu_C3';
    curTable.(name) =otsuAll(:,2).*256;
    name = 'std_C2';
    curTable.(name) =stdAll(:,1);
    name = 'std_C3';
    curTable.(name) =stdAll(:,2);
    name = 'good_C1_ID';
    curTable.(name) =passAll;
    
    %% cat table into final table
    if isempty(tableFinal)==1
        tableFinal=curTable;
    else
        tableFinal=[tableFinal;curTable];
    end

end
writetable(tableFinal,[tableName '_thresh.xlsx'])
