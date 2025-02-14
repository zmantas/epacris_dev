% check and compare Gibbs free energy of formation

clear all
close all

Gibbs_Chem = '../../Library/Gibbs/GibbsForm.dat';
g_ch = importdata(Gibbs_Chem);
temp = g_ch(1,:);

Gibbs_Eq = '../../Library/ChemEqu/molecules_all_forchecking.dat';
g_eq1 = importdata(Gibbs_Eq,'\t',3);
g_eq = g_eq1.data;

STD = 46;

for i=1:length(g_eq)
    if g_eq(i,1)==STD
        a=g_eq(i,2);
        b=g_eq(i,3);
        c=g_eq(i,4);
        d=g_eq(i,5);
        e=g_eq(i,6);
        C=g_eq(i,9);
        H=g_eq(i,7);
        N=g_eq(i,10);
        O=g_eq(i,11);
        S=g_eq(i,19);
    end
end
g = (a./temp+b+c*temp+d*temp.^2+e*temp.^3)*4.1840;
for i=1:length(g_eq)
    if g_eq(i,1)==53
        a=g_eq(i,2);
        b=g_eq(i,3);
        c=g_eq(i,4);
        d=g_eq(i,5);
        e=g_eq(i,6);
    end
end
g_H2 = (a./temp+b+c*temp+d*temp.^2+e*temp.^3)*4.1840;
for i=1:length(g_eq)
    if g_eq(i,1)==54
        a=g_eq(i,2);
        b=g_eq(i,3);
        c=g_eq(i,4);
        d=g_eq(i,5);
        e=g_eq(i,6);
    end
end
g_O2 = (a./temp+b+c*temp+d*temp.^2+e*temp.^3)*4.1840;
for i=1:length(g_eq)
    if g_eq(i,1)==55
        a=g_eq(i,2);
        b=g_eq(i,3);
        c=g_eq(i,4);
        d=g_eq(i,5);
        e=g_eq(i,6);
    end
end
g_N2 = (a./temp+b+c*temp+d*temp.^2+e*temp.^3)*4.1840;

g_c = g_ch(STD+1,:);
g_c = g_c - g_ch(19+1,:)*C - g_ch(40+1,:)*S;
plot(temp,g_c,'LineWidth',2);
hold on

g = g-0.5*g_H2*H-0.5*g_O2*O-0.5*g_N2*N;
plot(temp,g,'LineWidth',2);

%% correct the one that is incorrect
%% in this case it was 46 (HS)
% g_correct = g + g_ch(19+1,:)*C + g_ch(40+1,:)*S;