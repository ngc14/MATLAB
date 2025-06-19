function dist = dFunc(x,y)
x = x';
y = y';
dist = NaN(1,size(y,2));
for j = size(y,2):-1:1
    dist(j)= 1-xcorr(x,y(:,j),0,'normalized');
end
% dist = zeros(size(y,1),1);
% for i=1:size(y,1)
%     dist(i) = 1-corr(x',y(i,:)','Type','Spearman');
% end
end