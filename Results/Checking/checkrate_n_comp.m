% check reaction rates
clear all
close all

KB=1.38064852E-23;

%DIR = '../KEPLER167E_PhotoN/Planet_KEPLER167E_g25_meta2.0_i100/';
%DIR = '../PH2B_PhotoN/Planet_PH2B_g25_meta0.0_i100/';
DIR = '../WASP39B_PICASO_C1.0_R0.2_T100_K11_S/';
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

% Type
TYPE = 'R';
% STD
STD = [1,2,3,4,5,6,7,8];

% Reaction Rate
f=fopen([DIR,'ChemicalRate0.dat'],'r');
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

% Reaction Rate
f=fopen([DIR,'y/ChemicalRate0.dat'],'r');
delimiter = '\t';
startRow = 2;
formatSpec = '%s';
for j=1:zbin
        formatSpec = [formatSpec,'%f'];
end
formatSpec = [formatSpec,'%[^\n\r]'];
dataArray = textscan(f, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
ChemicalRate_comp = [dataArray{1:end-1}];
clearvars delimiter startRow formatSpec dataArray
fclose(f);

% Reaction
count=1;
for i=1:length(STD)
    for j=1:size(ChemicalRate,1)
            if ChemicalRate(j,1)==[TYPE,num2str(STD(i))]
                for k=1:size(ChemicalRate,2)-1
                    rate(k,count) = str2num(ChemicalRate{j,k+1});
                end
                Reaction{count}=ChemicalRate(j,1);
                count=count+1;
            end
    end
end

count=1;
for i=1:length(STD)
    for j=1:size(ChemicalRate,1)
            if ChemicalRate(j,1)==[TYPE,num2str(STD(i))]
                for k=1:size(ChemicalRate,2)-1
                    rate_comp(k,count) = str2num(ChemicalRate_comp{j,k+1});
                end
                %Reaction{count}=rname;
                count=count+1;
            end
    end
end


frate = sum(rate)*thickl;
[frate_s,idx]=sort(frate,'descend');

% Create figure
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
%for i=[1,3:length(Reaction)]
for i=idx
    loglog(abs(rate(:,i)),p,'DisplayName',[Reaction{i},' ',num2str(frate(i),'%.2e')],'LineWidth',2);
    loglog(abs(rate_comp(:,i)),p,'--','LineWidth',2);
end
%loglog(sum(max(0,rate),2),p,'k','DisplayName','Total Source','LineWidth',2);
%loglog(-sum(min(0,rate),2),p,'--k','DisplayName','Total Sink','LineWidth',2);
xlabel('Rate [cm-3 s-1]');
ylabel('Pressure [Pa]');
%xlim(axes1,[1E-6 1E+2]);
box(axes1,'on');
set(axes1,'FontSize',20,'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log','YDir','reverse');
legend(axes1,'show','location','BestOutside');
set(gcf,'Position',[0,0,900,500]);