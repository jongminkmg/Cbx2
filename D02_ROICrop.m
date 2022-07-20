% Function ROICrop
% : Given 1) a directory of images, 2) Coordinates [x, y, offset_direction]
% : Return cropped 13pixel x 13pixel rectangle images around the
% coordinates.
% : note, offset. Default, right(1). If outside of cell, then down->left->up.
% : 1) right  x+14
% : 2) down   y+14
% : 3) left   x-14
% : 4) up     y-14
% 
% You need to choose the directory first
% In the directory, you need 128x128 (or other sizes) cropped images.
% If one nuclei have more than one puncta, crop the nuclei again.
% Then get the x,y coordinates of the puncta center. 
% (Use coordinate picker tool in the ImageJ)
% The third column will be offset direction. 

% Last updated: 2022-07-20, Jongmin Kim (jongminkmg@gmail.com)
% Please note this function is intended for my own use, and has fixed input requirements. 
% Feel free to use/modify in any way.


function ROICrop (Crd)

Nuclei = dir ('*.tif');
n = length(Nuclei);
counter = 1;

for i=1:n
    nuc = imread(Nuclei(i).name);
    if (Crd(i,1) >= 6) && (Crd(i,2) >= 6) && (Crd(i,1)<=122) && (Crd(i,2)<=122) % chk if rectangle will not outside of the image
        nucCrop = imcrop (nuc, [Crd(i,1)-5 Crd(i,2)-5 12 12]);
        if (Crd(i,3) == 1) % offset value 1, move to right
            offCrop = imcrop (nuc, [Crd(i,1)-5+14 Crd(i,2)-5 12 12]);
        elseif (Crd(i,3) == 2) % offset value 2, move to down
            offCrop = imcrop (nuc, [Crd(i,1)-5 Crd(i,2)-5+14 12 12]);
        elseif (Crd(i,3) == 3) % offset value 3, move to left
            offCrop = imcrop (nuc, [Crd(i,1)-5-14 Crd(i,2)-5 12 12]);
        elseif (Crd(i,3) == 4) % offset value 3, move to left
            offCrop = imcrop (nuc, [Crd(i,1)-5 Crd(i,2)-5-14 12 12]);
        end
        %image(nuc)
        %image(nucCrop)
        imwrite(nucCrop, strcat('Puncta_C',int2str(counter), '_', Nuclei(i).name))
        imwrite(offCrop, strcat('Offset_C',int2str(counter), '_', Nuclei(i).name))
        
        if counter >= 4 % change this number based on the number of channels you have
           counter = 1;
        else
           counter = counter + 1;
        end
    end
end
