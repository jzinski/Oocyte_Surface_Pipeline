function [minRadBW_CH_Tot]=combineBinaries(minRadBW_CH, minRadBW_CHp,imstksm,individOocyteOutlines,baseName)

marker1=min(find(max(max(minRadBW_CH))==1));
marker2=max(find(max(max(minRadBW_CH))==1));
minRadBW_CH_Tot=zeros(size(minRadBW_CH,2),size(minRadBW_CH,1),size(minRadBW_CH,3));
minRadBW_CH_Tot(:,:,1:marker1)=minRadBW_CH(:,:,1:marker1) | minRadBW_CHp(:,:,1:marker1);
minRadBW_CH_Tot(:,:,marker1:marker2)=minRadBW_CH(:,:,marker1:marker2) ;
minRadBW_CH_Tot(:,:,marker2:size(minRadBW_CH,3))=minRadBW_CH(:,:,marker2:size(minRadBW_CH,3)) | minRadBW_CHp(:,:,marker2:size(minRadBW_CH,3));

CH2=uint8(zeros(size(minRadBW_CH_Tot,1),size(minRadBW_CH_Tot,2),3,size(minRadBW_CH_Tot,3)));
for i=1:numel(minRadBW_CH(1,1,:))
    CH2(:,:,:,i)=imfuse(imstksm(:,:,i),minRadBW_CH_Tot(:,:,i));
end
CH2=permute(CH2,[1,2,4,3]);
if individOocyteOutlines==1
figure
s=sliceViewer(CH2);
    
    hAx = getAxesHandle(s);

    filename = [baseName '_c1.gif'];

    sliceNums = 1:numel(minRadBW_CH(1,1,:));


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

end