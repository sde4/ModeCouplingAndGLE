modset='10x10/';
area = 'Lx0.04Ly0.04/';
%area = 'Lx1.00Ly1.00/';
f1=strcat('/Users/Subhadeep_De/Taub/taub10/AREND/CG/ModeCouplingAndGLE/', modset, area,'alphamatsorted.dat');
alphamat = importdata(f1);
f1=strcat('/Users/Subhadeep_De/Taub/taub10/AREND/CG/ModeCouplingAndGLE/', modset, area,'IRs_pqrcombmatsorted.dat');
IR = importdata(f1);
f1=strcat('/Users/Subhadeep_De/Taub/taub10/AREND/CG/ModeCouplingAndGLE/', modset, area,'modindmat.dat');
modindmat2=importdata(f1);
f1=strcat('/Users/Subhadeep_De/Taub/taub10/AREND/CG/ModeCouplingAndGLE/', modset, area,'frvec.dat');
om2=2*pi*importdata(f1);
f1=strcat('/Users/Subhadeep_De/Taub/taub10/AREND/CG/ModeCouplingAndGLE/', modset, area,'IRs_pqrcountmat.dat');
IR2countmat = importdata(f1);


[val, ind] = sortrows(IR(:, 1:4));
IR2 = IR(ind, :);
alphamat2 = alphamat(ind, :);
IR2(:, 5) = alphamat2(:, 1);
IR2(:, 6) = [];

t1=strtrim(cellstr(strcat(num2str(modindmat2(IR2(:, 2), 1)),',' ,num2str(modindmat2(IR2(:, 2), 2))))')';
t2=strtrim(cellstr(strcat(num2str(modindmat2(IR2(:, 2), 1)),',' ,num2str(modindmat2(IR2(:, 2), 2))))')';
t3=strtrim(cellstr(strcat(num2str(modindmat2(IR2(:, 3), 1)),',' ,num2str(modindmat2(IR2(:, 3), 2))))')';
t4=strtrim(cellstr(strcat(num2str(modindmat2(IR2(:, 4), 1)),',' ,num2str(modindmat2(IR2(:, 4), 2))))')';
t5=strtrim(cellstr(num2str(IR2(:, 5)))')';
IR2str=[t1, t2, t3, t4, t5];


IR2(:, 6) = om2(IR2(:, 1));
IR2(:, 7) = om2(IR2(:, 2));
IR2(:, 8) = om2(IR2(:, 3));
IR2(:, 9) = om2(IR2(:, 4));
IR2(:, 10) = IR2(:, 5)./(IR2(:, 6).*IR2(:, 7).*IR2(:, 8).*IR2(:, 9));

%figure; hold on;
%for ii=1:36
%    modeid=ii;plot(cumsum(IR2(IR2countmat(modeid, 3)-IR2countmat(modeid, 2)+1:IR2countmat(modeid, 3), 10)), '.')
%end
