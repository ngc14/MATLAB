%  Read ivf file (image formate used in Optical Imaging)
%  [image, width, height]=ReadIVF(filename)

function image=ReadIVF(filename)

%filename='c:\_\testivf.ivf';
fidimg=fopen(filename, 'r');
header=fread(fidimg, 3, 'long');
if header(1)~=4
    fprintf('Sorry, %s is not a valid IVF file (type~=4)!\n', filename);
end
width=header(2);
height=header(3);
image=fread(fidimg, width*height, 'float');
image=permute(reshape(image, [width, height]), [2,1]);
fclose(fidimg);

return;
