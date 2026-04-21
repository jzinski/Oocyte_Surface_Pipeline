function [imstksm]=smooth3D(imstk, x1, y1, z1)
%outputs a 3-D smoothed stack in identical format
%imstk=stack of incoming images  (x by y by z)
%x1=diamteter of cell in x direction  (usually 9)
%y1=diameter of cell in y direction (usually 9)
%z1=diameter of cell in z direction (usually 3)


imstk8=uint8(imstk);



%smooths%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%diskfilter:
%smoothdapi8=imfilter(dapi8, fspecial('disk', radius));

h = fspecial3_mod('gaussian', [x1,y1,z1],[x1,y1,z1]);   %generates a 3-D kernel

imstksm=imfilter(imstk,h);   %filers in 3-D


