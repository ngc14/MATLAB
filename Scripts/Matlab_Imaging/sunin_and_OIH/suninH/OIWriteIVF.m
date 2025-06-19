%  Write ivf file (image formate used in Optical Imaging
%  status=WriteIVF(data, filename)

function status=WriteIVF(data, filename)

[height, width]=size(data);
fid=fopen(filename, 'w');
image=flipud(imrotate(data,90));
image=reshape(image, 1, height*width);
header=[4, width, height];
fwrite(fid, header, 'long');
status=fwrite(fid, image, 'float');
fclose(fid);

return;
