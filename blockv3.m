% Blocking method for calculating standard deviation of 'correlated data'
% Reference: H. Flyvbjerg and H. G. Petersen, JCP, 91, 461 (1989).
% Jaeoh Shin, 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load eteFile.txt
load rg.txt
c1=eteFile(:,2); % two source data-contact number and overlap length
c111=rg(:,2);

ll0=length(c1);
trial=floor(log2(ll0))-1; % the number of transformation.


for i=1:trial  % For each loop, the data size becomes half of previous one.

ll(i)=length(c1);
var(i)=std(c1,1)*std(c1,1)/(ll(i)-1); % variance of the data- Eq. (26) of the reference
var111(i)=std(c111,1)*std(c111,1)/(ll(i)-1); % variance of the data- Eq. (26) of the reference
cc=[];
cc111=[];
M=mod(ll(i), 2); % 0 or 1
for k=1:2:ll(i)-M
    cc=[cc (c1(k)+c1(k+1))/2.0];
    cc111=[cc111 (c111(k)+c111(k+1))/2.0];
end
c1=cc;
c111=cc111;
end



% Writhe the standard deviation and its errorbar for each trials
fname=sprintf('%s', 'error.dat');
fout=fopen(fname,'wt');

fname111=sprintf('%s', 'error111.dat');
fout111=fopen(fname111,'wt');

for k=1:trial
    fprintf(fout, '%d\t%f\t%f\n', k, sqrt(var(k)), sqrt(var(k)/(2*ll(k)-1)))
    fprintf(fout111, '%d\t%f\t%f\n', k, sqrt(var111(k)), sqrt(var111(k)/(2*ll(k)-1)))
end
