function [I] = FRavg(F)
% F:the spectrum of an image
% I:the average intensity in each R
F=abs(F);
[d1,d2]=size(F);
x=1:d1;
y=1:d2;
x=x-fix(d1/2);
y=y-fix(d2/2);
[X,Y]=meshgrid(x,y);
R=sqrt(X.^2+Y.^2);
R=fix(R);
saver=zeros(1,d1/2);
% Seeking average
for i=0:fix(d1/2)-1
index=find(R==i);
[row, col] = ind2sub(size(R), index);
count=0;
for j=1:length(row)
    count=count+F(row(j),col(j));
end
count=count/length(row);
saver(i+1)=count;
end
I=saver;
end

