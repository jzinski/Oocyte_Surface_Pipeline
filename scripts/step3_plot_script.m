%% Path setup - adds ../functions and ../thirdparty to MATLAB path
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir,'..','functions'));
addpath(genpath(fullfile(scriptDir,'..','thirdparty')));

%% Data folder
dataDir = fullfile(scriptDir,'..','data');
cd(dataDir);

%% user input
tableName='Oocyte_sample_list_thresh2';

%% x y z pixel size microns (user can change)
xsize=0.3321;
ysize=0.3321;
zsize=2.6;

qualThresh=0;

%% color for each genotype
genoColorSet=[.75,0,0;0,.75,0;0,0,.75];
dateColorSet=[.75,0.75,0.75;0.5,.5,0.5;0.25,0.25,.25;0.75,0,.75];
%% get thresholds for C2 and c3
metaAll=readtable([tableName '.xlsx']);
metaAllPrune=metaAll(metaAll.good_C1_ID>qualThresh & metaAll.C3_Thresh~=0,:);

%% compile variables 

for i=1:numel(metaAllPrune.name)
    load([metaAllPrune.name{i} '_out.mat'])
    elePercC2_all(i,:)=elePercC2;
    elePercC3_all(i,:)=elePercC3;
    H1counts_all(i,:)=H1counts;
    H1countsC2_all(i,:)=H1countsC2;
    Cencounts_all(i,:)=Cencounts;
    i
end
H1volume_all=(H1counts_all)./Cencounts_all;  %normalizing for total amount of pixels
H1volumeC2_all=(H1countsC2_all)./Cencounts_all;  %normalizing for total amount of pixels

%% Plot C2 percentage Coverage
uniqueDates=unique(metaAllPrune.date);
uniqueGeno=unique(metaAllPrune.genotype);

for i=1:numel(uniqueDates)

    curDateIdx=metaAllPrune.date==uniqueDates(i);
    figure
    try
        subplot(2,2,1)
        hold on
        plot(0:5:175,elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:),'-','LineWidth', 1, 'Color', [.75,.25,.25])  
        plot(0:5:175,elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:),'-','LineWidth', 1, 'Color', [.25,.75,.25])   
        plot(0:5:175,elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:),'-','LineWidth', 1, 'Color', [.25,.25,.75])  
        h1=plot(0:5:175,mean(elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:)),'LineWidth', 4, 'Color', [.75,0,0]);  
        h2=plot(0:5:175,mean(elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:)),'LineWidth', 4, 'Color', [0,.75,0]); 
        h3=plot(0:5:175,mean(elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:)),'LineWidth', 4, 'Color', [0,0,.75]);
        
        shadedErrorBar(0:5:175,mean(elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:)), std(elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:)),'lineProps','b')

        legend([h1,h2,h3],uniqueGeno)
        xlim([0,180])
        ylim([0,1])
        xlabel('Angle From Animal Pole (degrees)', 'FontSize', 12);
        ylabel('Fraction Surface Covered by Dazl', 'FontSize', 12);
        title(['CyclinB1 Fraction Coverage'])
        

        curDateIdx=metaAllPrune.date==uniqueDates(i);
        subplot(2,2,2)
        hold on
        plot(0:5:175,elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:),'LineWidth', 1, 'Color', [.75,0,0])  
        plot(0:5:175,elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:),'LineWidth', 1, 'Color', [0,.75,0])   
        plot(0:5:175,elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:),'LineWidth', 1, 'Color', [0,0,.75])   
        plot(0:5:175,mean(elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:)),'LineWidth', 4, 'Color', [.75,0,0])  
        plot(0:5:175,mean(elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:)),'LineWidth', 4, 'Color', [0,.75,0])   
        plot(0:5:175,mean(elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:)),'LineWidth', 4, 'Color', [0,0,.75])  
        xlim([0,180])
        xlabel('Angle From Animal Pole (degrees)', 'FontSize', 12);
        ylabel('Fraction Surface Covered by Dazl', 'FontSize', 12);
        title(['Dazl Fraction Coverage'])


        subplot(2,2,3)
        hold on
        plot(0:5:175,H1volumeC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:),'LineWidth', 1, 'Color', [.75,0,0])  
        plot(0:5:175,H1volumeC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:),'LineWidth', 1, 'Color', [0,.75,0])   
        plot(0:5:175,H1volumeC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:),'LineWidth', 1, 'Color', [0,0,.75])   
        plot(0:5:175,mean(H1volumeC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:)),'LineWidth', 4, 'Color', [.75,0,0])  
        plot(0:5:175,mean(H1volumeC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:)),'LineWidth', 4, 'Color', [0,.75,0])   
        plot(0:5:175,mean(H1volumeC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:)),'LineWidth', 4, 'Color', [0,0,.75])  
        grid on;
        xlim([0, 180]);
        xlabel('Angle From Animal Pole (degrees)', 'FontSize', 12);
        ylabel('Percent Volume C2', 'FontSize', 12);

        subplot(2,2,4)
        hold on
        plot(0:5:175,H1volume_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:),'LineWidth', 1, 'Color', [.75,0,0])  
        plot(0:5:175,H1volume_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:),'LineWidth', 1, 'Color', [0,.75,0])   
        plot(0:5:175,H1volume_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:),'LineWidth', 1, 'Color', [0,0,.75])   
        plot(0:5:175,mean(H1volume_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:)),'LineWidth', 4, 'Color', [.75,0,0])  
        plot(0:5:175,mean(H1volume_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:)),'LineWidth', 4, 'Color', [0,.75,0])   
        plot(0:5:175,mean(H1volume_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:)),'LineWidth', 4, 'Color', [0,0,.75])  
        grid on;
        xlim([0, 180]);
        xlabel('Angle From Animal Pole (degrees)', 'FontSize', 12);
        ylabel('Percent Volume C3', 'FontSize', 12);

        sgtitle(['Date = ' num2str(uniqueDates(i))])
        set(gcf,'units','points','position',[10,10,900,900])
    catch
        sgtitle(['Date = ' num2str(uniqueDates(i))])
        set(gcf,'units','points','position',[10,10,900,900])
    end
end


%% Plot first two dates for pub
uniqueDates=unique(metaAllPrune.date);
uniqueGeno=unique(metaAllPrune.genotype);

   curDateIdx=metaAllPrune.date==uniqueDates(1) | metaAllPrune.date==uniqueDates(2);
    figure

        hold on
     %   plot(0:5:175,elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:),'-','LineWidth', 1, 'Color', [1,.75,.75])  
    %    plot(0:5:175,elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:),'-','LineWidth', 1, 'Color', [.75,1,.75])   
   %     plot(0:5:175,elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:),'-','LineWidth', 1, 'Color', [.75,.75,1])  

        shadedErrorBar(2.5:5:177.5,mean(elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:)), std(elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:)),'lineProps','r')
        shadedErrorBar(2.5:5:177.5,mean(elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:)), std(elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:)),'lineProps','g')
        shadedErrorBar(2.5:5:177.5,mean(elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:)), std(elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:)),'lineProps','b')
        plot(2.5:5:177.5,mean(elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:)), 'LineWidth', 3, 'Color', [.75,0,0])
        plot(2.5:5:177.5,mean(elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:)), 'LineWidth', 3, 'Color', [0,.75,0])
        plot(2.5:5:177.5,mean(elePercC2_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:)), 'LineWidth', 3, 'Color', [0,0,.75])

        
        set(gca,'FontSize',12)
        set(gca,'LineWidth',2)
        xlim([0,180])
        ylim([0,1])
        xlabel('Angle From Animal Pole (degrees)', 'FontSize', 14);
        ylabel('Fraction Surface Covered by CyclinB1', 'FontSize', 14);
        title(['CyclinB1 Fraction Coverage'],'FontSize', 14)
        grid on
        set(gcf,'units','points','position',[100,100,500,300])

        figure
        hold on
     %   plot(0:5:175,elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:),'LineWidth', 1, 'Color', [1,.75,.75])  
      %  plot(0:5:175,elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:),'LineWidth', 1, 'Color', [.75,1,.75])   
       % plot(0:5:175,elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:),'LineWidth', 1, 'Color', [.75,.75,1])   
 
        
        
        shadedErrorBar(2.5:5:177.5,mean(elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:)), std(elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:)),'lineProps','r')
        shadedErrorBar(2.5:5:177.5,mean(elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:)), std(elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:)),'lineProps','g')
        shadedErrorBar(2.5:5:177.5,mean(elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:)), std(elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:)),'lineProps','b')

        plot(2.5:5:177.5,mean(elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(1)),:)), 'LineWidth', 3, 'Color', [.75,0,0])
        plot(2.5:5:177.5,mean(elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(2)),:)), 'LineWidth', 3, 'Color', [0,.75,0])
        plot(2.5:5:177.5,mean(elePercC3_all(curDateIdx & strcmp(metaAllPrune.genotype,uniqueGeno(3)),:)), 'LineWidth', 3, 'Color', [0,0,.75])

        set(gca,'FontSize',12)
        set(gca,'LineWidth',2)
        xlim([0,180])
        ylim([0,1])
        xlabel('Angle From Animal Pole (degrees)', 'FontSize', 14);
        ylabel('Fraction Surface Covered by Dazl', 'FontSize', 14);
        title(['Dazl Fraction Coverage'],'FontSize', 14)
        grid on



        set(gcf,'units','points','position',[600,100,500,300])


%% plot all dates on same plots

figure
for i=1:numel(uniqueGeno)
   
       curGenoIdx=strcmp(metaAllPrune.genotype,uniqueGeno(i));
        subplot(3,2,i*2-1)
        hold on
        for o=1:numel(uniqueDates)
           curDateIdx=metaAllPrune.date==uniqueDates(o);
           plot(0:5:175,mean(elePercC2_all(curDateIdx & curGenoIdx,:)),'LineWidth', 3, 'Color', dateColorSet(o,:) );          
        end
        if i==1
            legend(cellstr(num2str(uniqueDates)))
        end
        xlim([0,180])
        xlabel('Angle From Animal Pole (degrees)', 'FontSize', 12);
        ylabel('Fraction Surface Covered by Dazl', 'FontSize', 12);
        title(['CyclinB1 Fraction Coverage ' uniqueGeno(i)])

        subplot(3,2,i*2)
        hold on
        for o=1:numel(uniqueDates)
           curDateIdx=metaAllPrune.date==uniqueDates(o);
           plot(0:5:175,mean(elePercC3_all(curDateIdx & curGenoIdx,:)),'LineWidth', 3, 'Color', dateColorSet(o,:) );          
        end
        xlim([0,180])
        xlabel('Angle From Animal Pole (degrees)', 'FontSize', 12);
        ylabel('Fraction Surface Covered by Dazl', 'FontSize', 12);
        title(['Dazl Fraction Coverage ' uniqueGeno(i)])
        
        sgtitle(['Date Compare by genotype'])
        set(gcf,'units','points','position',[10,10,900,900])
end

%% fit sigmoid to each date and geno
C2bounds_output=zeros(numel(metaAllPrune.name),1);
C3bounds_output=zeros(numel(metaAllPrune.name),1);
ft = fittype('c3/(1 + exp(-c1*(x-c2)))','indep','x');
ftn = fittype('c3/(1 + exp(-c1*(-x+c2)))','indep','x');
clear C2boundAll C3boundAll C2boundsAllc1 C3boundsCurc1
count=0;
figure
hold on
for i=1:numel(uniqueGeno)
    curGenoIdx=strcmp(metaAllPrune.genotype,uniqueGeno(i));
    for o=1:numel(uniqueDates)
       curDateIdx=metaAllPrune.date==uniqueDates(o);
       curFitPtsC2=elePercC2_all(curDateIdx & curGenoIdx,:);     
       curFitPtsC3=elePercC3_all(curDateIdx & curGenoIdx,:);
       count=count+1;
       clear C2boundsCur C2boundsCurc1 C3boundsCur C3boundsCurc1
        subplot(numel(uniqueGeno),numel(uniqueDates),count)
        hold on
       for u=1:sum(curDateIdx & curGenoIdx)
           C2fitCur=fit([0:5:175]',curFitPtsC2(u,:)',ft, 'StartPoint', [-.5,30,1], 'Lower',[-1,5,0],'Upper',[-.09,70,1]);
           plot(0:5:175,curFitPtsC2(u,:),'LineWidth', 3, 'Color', dateColorSet(o,:) );    
           plot(C2fitCur)
           C2boundsCur(u)=C2fitCur.c2;
           C2boundsCurc1(u)=C2fitCur.c1;
           C2boundAll{i,o}=C2boundsCur;
           C2boundsAllc1{i,o}=C2boundsCurc1;

           
           C3fitCur=fit([0:5:175]',curFitPtsC3(u,:)',ftn, 'StartPoint', [-.2,120,1]);
           plot(0:5:175,curFitPtsC3(u,:),'LineWidth', 3, 'Color', dateColorSet(o,:) );    
           plot(C3fitCur)
           C3boundsCur(u)=C3fitCur.c2;
           C3boundsCurc1(u)=C3fitCur.c1;
           C3boundAll{i,o}=C3boundsCur;   
           C3boundsAllc1{i,o}=C3boundsCurc1;
           legend('off')
       
       end
       if isempty(curFitPtsC2)==0
       C2bounds_output(curDateIdx & curGenoIdx)=C2boundsCur;
       C3bounds_output(curDateIdx & curGenoIdx)=C3boundsCur;
       end
    end
end


table(metaAllPrune.genotype,metaAllPrune.date,C2bounds_output,C3bounds_output,'VariableNames', {'1','2','3','4'});
C2bounds_output_full=zeros(size(metaAll,1),1);
C2bounds_output_full(metaAll.good_C1_ID>qualThresh & metaAll.C3_Thresh~=0)=C2bounds_output;
C3bounds_output_full=zeros(size(metaAll,1),1);
C3bounds_output_full(metaAll.good_C1_ID>qualThresh & metaAll.C3_Thresh~=0)=C3bounds_output;

elePercC2_all_full=zeros(size(metaAll,1),size(elePercC2_all,2));
elePercC2_all_full(metaAll.good_C1_ID>qualThresh & metaAll.C3_Thresh~=0,:)=elePercC2_all;


elePercC3_all_full=zeros(size(metaAll,1),size(elePercC3_all,2));
elePercC3_all_full(metaAll.good_C1_ID>qualThresh & metaAll.C3_Thresh~=0,:)=elePercC3_all;


% makes a tibble table of angles for C2 and C3 ercentages
angleCount=[0:5:175]';
for i=1:size(elePercC3_all_full,1)
    curElesC2=elePercC2_all_full(i,:)';    
    curElesC3=elePercC3_all_full(i,:)';  
    qual_Cur=metaAll.good_C1_ID(i)*ones(size(curEles)); 
    
    date_Cur=metaAll.date(i)*ones(size(curEles)); 
    genoCur=repmat(metaAll.genotype(i),size(curEles));   
    nameCur=repmat(metaAll.name(i),size(curEles));
    
    if i==1
        eleTable=table(nameCur,genoCur,date_Cur,qual_Cur,angleCount,curElesC2,curElesC3,'VariableNames', {'Name','Condition','date','Quality','Angle_From_Animal','ElePercC2','ElePercC3'});
    else
        eleTableNew=table(nameCur,genoCur,date_Cur,qual_Cur,angleCount,curElesC2,curElesC3,'VariableNames', {'Name','Condition','date','Quality','Angle_From_Animal','ElePercC2','ElePercC3'});
        eleTable=cat(1,eleTable,eleTableNew);
    end
    
end

writetable(eleTable,[tableName '_elePerc.xlsx'])

metaAll.C2_Bound=C2bounds_output_full;
metaAll.C3_Bound=C3bounds_output_full;

writetable(metaAll,[tableName '_2.xlsx'])
writematrix(elePercC2_all_full,[tableName 'elePercC2.xlsx'])
writematrix(elePercC3_all_full,[tableName 'elePercC3.xlsx'])

%% bounds and volumes
count=0;
figure 
hold on
for o=1:size(C2boundAll,2)
    for i=1:size(C2boundAll,1)
        count=count+1;
        plot((rand(size(C2boundAll{i,o}))*.2)+ones(size(C2boundAll{i,o}))*count-.4, C2boundAll{i,o},'.','Color',[genoColorSet(i,:)])
     %  violinplot(C2boundAll{i,o})
        ylim([0,180])
        xlim([0,numel(C2boundAll)+1])
    end
end

count=0;
figure 
hold on
for o=1:size(C3boundAll,2)
    for i=1:size(C3boundAll,1)
        count=count+1;
        plot((rand(size(C3boundAll{i,o}))*.2)+ones(size(C3boundAll{i,o}))*count-.4, C3boundAll{i,o},'.','Color',[genoColorSet(i,:)])
     %  violinplot(C2boundAll{i,o})
        ylim([0,180])
        xlim([0,numel(C3boundAll)+1])
    end
end

figure 
hold on
for o=1:2
    for i=1:size(C3boundAll,1)
        count=count+1;
        plot(C2boundAll{i,o}, C3boundAll{i,o},'.','Color',[genoColorSet(i,:)])
        hold on
     %  violinplot(C2boundAll{i,o})
        ylim([0,180])
        xlim([0,180])
    end
end

figure 
hold on
for o=1:2
    for i=1:size(C3boundAll,1)
        curDateIdx=metaAllPrune.date==uniqueDates(o);
        curGenoIdx=strcmp(metaAllPrune.genotype,uniqueGeno(i));

        count=count+1;
        plot(metaAllPrune.Volume_C1(curGenoIdx & curDateIdx), C3boundAll{i,o},'.','Color',[genoColorSet(i,:)])
        hold on
     %  violinplot(C2boundAll{i,o})
        ylim([0,180])
      % xlim([0,180])
    end
end

figure 
hold on
for o=1:2
    for i=1:size(C3boundAll,1)
        curDateIdx=metaAllPrune.date==uniqueDates(o);
        curGenoIdx=strcmp(metaAllPrune.genotype,uniqueGeno(i));

        count=count+1;
        plot(metaAllPrune.Volume_C1(curGenoIdx & curDateIdx).^(1/3), C3boundAll{i,o},'.','Color',[genoColorSet(i,:)])
        hold on
     %  violinplot(C2boundAll{i,o})
        ylim([0,180])
      % xlim([0,180])
    end
end

figure 
hold on
for o=1:2
    for i=1:size(C3boundAll,1)
        curDateIdx=metaAllPrune.date==uniqueDates(o);
        curGenoIdx=strcmp(metaAllPrune.genotype,uniqueGeno(i));

        count=count+1;
        plot(metaAllPrune.Volume_C2(curGenoIdx & curDateIdx), C2boundAll{i,o},'.','Color',[genoColorSet(i,:)])
        hold on
     %  violinplot(C2boundAll{i,o})
        ylim([0,180])
      % xlim([0,180])
    end
end

figure 
hold on
for o=1:2
    for i=1:size(C3boundAll,1)
        curDateIdx=metaAllPrune.date==uniqueDates(o);
        curGenoIdx=strcmp(metaAllPrune.genotype,uniqueGeno(i));

        count=count+1;
        plot(metaAllPrune.Volume_C3(curGenoIdx & curDateIdx), C3boundAll{i,o},'.','Color',[genoColorSet(i,:)])
        hold on
     %  violinplot(C2boundAll{i,o})
        ylim([0,180])
      % xlim([0,180])
    end
end

figure 
hold on
for o=1:2
    for i=1:size(C3boundAll,1)
        curDateIdx=metaAllPrune.date==uniqueDates(o);
        curGenoIdx=strcmp(metaAllPrune.genotype,uniqueGeno(i));

        count=count+1;
        plot(metaAllPrune.Volume_C3(curGenoIdx & curDateIdx), C3boundAll{i,o},'.','Color',[genoColorSet(i,:)])
        hold on
     %  violinplot(C2boundAll{i,o})
        ylim([0,180])
      % xlim([0,180])
    end
end

figure 
hold on
for o=1:2
    for i=1:size(C3boundAll,1)
        curDateIdx=metaAllPrune.date==uniqueDates(o);
        curGenoIdx=strcmp(metaAllPrune.genotype,uniqueGeno(i));

        count=count+1;
        plot(metaAllPrune.Volume_C3(curGenoIdx & curDateIdx),metaAllPrune.Volume_C3(curGenoIdx & curDateIdx), C3boundAll{i,o},'.','Color',[genoColorSet(i,:)])
        hold on
     %  violinplot(C2boundAll{i,o})
        ylim([0,180])
      % xlim([0,180])
    end
end

%% fit sigmoid to each date and geno
figure
hold on
count=0;
for i=1:numel(uniqueGeno)
    curGenoIdx=strcmp(metaAllPrune.genotype,uniqueGeno(i));
    for o=1:numel(uniqueDates)
       curDateIdx=metaAllPrune.date==uniqueDates(o);
       curFitPtsC2=elePercC2_all(curDateIdx & curGenoIdx,:);     
       curFitPtsC3=elePercC3_all(curDateIdx & curGenoIdx,:);
              count=count+1;
       clear C2boundsCur C3boundsCur
        subplot(numel(uniqueGeno),numel(uniqueDates),count)
        hold on
       for u=1:sum(curDateIdx & curGenoIdx)
             plot(0:5:175,curFitPtsC2(u,:),'LineWidth', 1, 'Color', [0,1,0]/2 );    
               
            plot(0:5:175,curFitPtsC3(u,:),'LineWidth', 1, 'Color', [1,0,0]./2 );    
              legend('off')
         xlim([0,180])
        xlabel('Angle From Animal Pole (degrees)', 'FontSize', 12);
        ylabel('Fraction Surface Covered by Dazl', 'FontSize', 12);
        title(['Fraction Coverage ' uniqueDates(o) ' ' uniqueGeno(i)])
        

       end
        plot(0:5:175,mean(curFitPtsC2),'LineWidth', 3, 'Color', [0,1,0] );    

        plot(0:5:175,mean(curFitPtsC3),'LineWidth', 3, 'Color', [1,0,0] );
    end
end
        set(gcf,'units','points','position',[10,10,900,900])
        

%% put extents into table

metaOutCat=metaAll;  %will hold new meta data generated by C2 C3 thresh function
blankAddT=table(zeros(size(metaAll,1),1),zeros(size(metaAll,1),1),zeros(size(metaAll,1),1),zeros(size(metaAll,1),1),zeros(size(metaAll,1),1),zeros(size(metaAll,1),1),zeros(size(metaAll,1),1),'VariableNames', {'Volume_C1','Volume_C2','Volume_C3','Mean_Int_C2','Mean_Int_C3', 'C2_Thresh', 'C3_Thresh'});
metaOutCat=[metaOutCat,blankAddT];  %initilize columns for new variables;