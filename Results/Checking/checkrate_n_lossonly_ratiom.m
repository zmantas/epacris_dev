% check reaction rates
clear all
close all

KB=1.38064852E-23;

%DIR = '../KEPLER167E_PhotoB/Planet_KEPLER167E_g25_meta1.0_i100/';
%DIR = '../PH2B_PhotoN/Planet_PH2B_g25_meta2.0_i100_deep/';
%DIR = '../K2_18B_Photo/Planet_K2_18B_g12_meta0.0_i60/';
%DIR = '../K2_18B_PhotoN/Planet_K2_18B_g12_meta0.0_i60/';
%DIR = '../WASP39B_Photo_NOS/Planet_WASP39B_g4_meta0.0_i100_k8_s_test/';
DIR = '../WASP39B_ELSIE_C1_morning/';

IM=importdata([DIR,'ConcentrationSTD.dat'],'\t',2);
B1=IM.data;
z = B1(1:end,1);
p = B1(1:end,5);
t = B1(1:end,4);
zbin = length(z);
thickl=(z(2)-z(1))*1e+5;
N1 = p./t/KB*1e-6;

% STD
STD = 47;

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
    if RP(i,1)==STD
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
    if RP(i,1)==STD || RP(i,2)==STD
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
    if RP(i,1)==STD || RP(i,2)==STD
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
    if RP(i,1)==STD
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

for i=1:length(Reaction)
    rate1(:,i)=rate(:,i)./B1(:,5+STD);
end

% compare with Kzz/H2
kzz = importdata([DIR,'KZZ.dat']);
kzz = 0.5*kzz(1:end-1,2)+0.5*kzz(2:end,2);
h = 1.38e-23*t/2.3/1.67e-27*1e+2;
tdiff = kzz./h.^2;

% Create figure
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
%for i=[1,3:length(Reaction)]
% for i=1:length(Reaction)
%     loglog(rate1(:,i),p,'DisplayName',[Reaction{i},' ',num2str(frate(i),'%.2e')],'LineWidth',2);
% end
loglog(sum(rate1(:,:),2),p,'DisplayName','Total Loss','LineWidth',4);
loglog(tdiff,p,'DisplayName','Diffusion','LineWidth',4);
xlabel('Rate [s-1]');
ylabel('Pressure [Pa]');
%xlim(axes1,[1E-6 1E+2]);
box(axes1,'on');
set(axes1,'FontSize',20,'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log','YDir','reverse');
legend(axes1,'show','location','BestOutside');
set(gcf,'Position',[0,0,900,500]);