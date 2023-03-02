clear all
load eteFile.txt


fig1=figure(1)
x=eteFile(:,2);
[c1,x1]=hist(x,160);
plot(x1,-log(c1)+4.7,'-');
xlabel ('End-to-end distance')
ylabel ('Free energy, F(x)/k_{B}T')

% fig2=figure(2)
% y=contact(:,1)/3600;
% [c2,x2]=hist(y,20);
% plot(x2,-log(c2),'r:o');
% xlabel ('Two chain overlap fraction, p')
% ylabel ('Free energy of overlap, F(p)/k_{B}T')