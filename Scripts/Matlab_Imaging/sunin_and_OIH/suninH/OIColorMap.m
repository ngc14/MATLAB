function colormap = OIColorMap(imin, clt)

% colormap = OIColorMap(imin, clt)
% input: imin is a uint8 map, clt is 256x3 matrix
% output: colormap is color map using clt

imin=double(imin)+1;
colormapt=zeros(size(imin, 1), size(imin,2), 3);
for i=1:size(imin,1)
    for j=1:size(imin, 2)
        colormapt(i,j,:)=clt(floor(imin(i,j)),:);
    end
end
colormap=norm_to_uint8(colormapt);
%colormap(:,:,1)=colormapt(:,:,1)*256;
%colormap(:,:,2)=colormapt(:,:,2)*256;
%colormap(:,:,3)=colormapt(:,:,3)*256;
colormap=uint8(colormap);
return