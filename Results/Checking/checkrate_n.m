% check reaction rates
clear all
close all

KB=1.38064852E-23;

%DIR = '../KEPLER167E_PhotoN/Planet_KEPLER167E_g25_meta2.0_i100/';
%DIR = '../PH2B_PhotoN/Planet_PH2B_g25_meta0.0_i100/';
DIR = '../WASP39B_ELSIE_C1_morning/';
%DIR = '../K2_18B_Photo_NOS/Planet_K2_18B_g12_meta0.0_i60_p5_lowCO2_noEmi/';
%DIR = '../JupiterN/Planet_JUPITER_nos_meta0.5_i100_strat_k3/';

IM=importdata([DIR,'ConcentrationSTD.dat'],'\t',2);
B1=IM.data;
z = B1(1:end,1);
p = B1(1:end,5);
t = B1(1:end,4);
zbin = length(z);
thickl=(z(2)-z(1))*1e+5;
N1 = p./t/KB*1e-6;

% STD
STD = 3;

% Reaction Rate
f=fopen([DIR,'ChemicalRate.dat'],'r');
delimiter = '\t';
startRow = 2;
formatSpec = '%s';
for j=1:zbin
        formatSpec = [formatSpec,'%f'];
end
formatSpec = [formatSpec,'%[^\n\r]'];
dataArray = textscan(f, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
ChemicalRate = [dataArray{1:end-1}];
clearvars delimiter startRow formatSpec dataArray
fclose(f);

% Reaction
count = 1;
% P
RP=importdata('../../Library/ReactionList/Reaction_P.txt');
for i=1:length(RP)
    if RP(i,1)==STD || RP(i,2)==STD || RP(i,3)==STD || RP(i,4)==STD
        for j=1:size(ChemicalRate,1)
            if ChemicalRate(j,1)==['P',num2str(i)]
                rname = ['P',num2str(i)];
                disp([rname,' ',num2str(RP(i,:))]);
                for k=1:size(ChemicalRate,2)-1
                    rate(k,count) = str2num(ChemicalRate{j,k+1});
                end
                Reaction{count}=rname;
                count=count+1;
            end
        end
    end
end
% R
RP=importdata('../../Library/ReactionList/Reaction_R.txt');
for i=1:length(RP)
    if RP(i,1)==STD || RP(i,2)==STD || RP(i,4)==STD || RP(i,5)==STD || RP(i,6)==STD
        for j=1:size(ChemicalRate,1)
            if ChemicalRate(j,1)==['R',num2str(i)]
                rname = ['R',num2str(i)];
                disp([rname,' ',num2str(RP(i,:))]);
                for k=1:size(ChemicalRate,2)-1
                    rate(k,count) = str2num(ChemicalRate{j,k+1});
                end
                Reaction{count}=rname;
                count=count+1;
            end
        end
    end
end
% M
RP=importdata('../../Library/ReactionList/Reaction_M.txt');
for i=1:length(RP)
    if RP(i,1)==STD || RP(i,2)==STD || RP(i,3)==STD || RP(i,4)==STD
        for j=1:size(ChemicalRate,1)
            if ChemicalRate(j,1)==['M',num2str(i)]
                rname = ['M',num2str(i)];
                disp([rname,' ',num2str(RP(i,:))]);
                for k=1:size(ChemicalRate,2)-1
                    rate(k,count) = str2num(ChemicalRate{j,k+1});
                end
                Reaction{count}=rname;
                count=count+1;
            end
        end
    end
end
% T
RP=importdata('../../Library/ReactionList/Reaction_T.txt');
for i=1:length(RP)
    if RP(i,1)==STD || RP(i,2)==STD || RP(i,3)==STD
        for j=1:size(ChemicalRate,1)
            if ChemicalRate(j,1)==['T',num2str(i)]
                rname = ['T',num2str(i)];
                disp([rname,' ',num2str(RP(i,:))]);
                for k=1:size(ChemicalRate,2)-1
                    rate(k,count) = str2num(ChemicalRate{j,k+1});
                end
                Reaction{count}=rname;
                count=count+1;
            end
        end
    end
end

%rate(:,1) = rate(:,1)+rate(:,2); 
frate = sum(rate)*thickl;
[frate_s,idx]=sort(frate,'descend');

% Create figure
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
%for i=[1,3:length(Reaction)]
for i=idx
    loglog(abs(rate(:,i)),p,'DisplayName',[Reaction{i},' ',num2str(frate(i),'%.2e')],'LineWidth',2);
end
%loglog(sum(max(0,rate),2),p,'k','DisplayName','Total Source','LineWidth',2);
%loglog(-sum(min(0,rate),2),p,'--k','DisplayName','Total Sink','LineWidth',2);
xlabel('Rate [cm-3 s-1]');
ylabel('Pressure [Pa]');
xlim(axes1,[1E-6 1E+2]);
box(axes1,'on');
set(axes1,'FontSize',20,'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log','YDir','reverse');
legend(axes1,'show','location','BestOutside');
set(gcf,'Position',[0,0,900,500]);