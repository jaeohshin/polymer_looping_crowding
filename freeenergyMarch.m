clear all
temp1=load('./eteFile.txt')

x1=temp1(:,2);
binnumber=100;
[c1,x11]=hist(x1,binnumber);

ll=length(c1)
norm=length(x1);
normi=1.0/norm;

fname=sprintf('%s','pdfuniform.txt');
fout=fopen(fname,'wt');

for k=1:ll
fprintf(fout,'{%f, %f}, ',x11(k), c1(k)*normi);
%fprintf(fout,'%f\t%e\n',x11(k), normi);
end
fprintf(fout,'\n\n');


