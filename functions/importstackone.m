function [stack] = importstackone(stackname, totchan, deschan)
%designed to input images from a single stack of tif images with 2
%channels.  it starts at im 1 and goes through the whole imageset

%stackname= the name of the image sequence up to the image numbers.  the
%for example: 'name_z'
%startim=number of first image in sequence. example: 1
%endim=number of last image in sequence. example: 20
%suffix=the channel suffix of desired images. example: '_c0002'

%stack= output of all images in x by y by z matrix

endim=numel(imfinfo([stackname '.tif']))/totchan;
startim=1;

startim=deschan+startim*totchan-totchan;
endim=deschan+endim*totchan-totchan;
j=1;
for i=startim:totchan:endim
    
    

stack(:,:,j)=imread([stackname '.tif'], 'index', i);


    j=j+1;
    
end

end