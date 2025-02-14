% Make Figure

%clear all;
%close all

% FILENAME = '../WASP39B_PICASO_C1.0_R0.2_T100_K11_S/ConcentrationSTD.dat';
FILENAME = '../k218b_soot/ConcentrationSTD_T.dat';

kb = 1.380658E-23;

% import data

IM=importdata(FILENAME,'\t',2);
B=IM.data;
N1 = sum(B(:,6:end)')';
z = B(1:end,1);
t = B(1:end,4);
p = B(1:end,5);
N = p./t/kb*1.0E-6;

B(:,end+1) = N-N1;

% select species

NAME = {'H2','H','He','H2O','OH','CH4','C2H2','CO','CO2','O2','O','O3','N2','NH3','NO','NO2','HCN','HCNO','S','S2','H2S','HS','SO','SO2','OCS','CS'};
		 
STD = [53,3,112,7,4,21,27,20,52,54,1,2,55,9,12,13,37,72,40,41,45,46,42,43,49,50];

SIZE = [6,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3];
		 
COLOR = {'k','--k',':k','k','-.k','g','-.g','--g',':g','b','--b',':b','m','--m',':m','-.m','c','--c','y','--y','r','--r',':r','-.r',':c','-.c'};

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'XScale','log','YScale','log','YDir','reverse',...
    'XTick',[1e-14 1e-12 1e-10 1e-08 1e-06 1e-4 1e-2 1],...    
    'XMinorTick','on','YMinorTick','on', 'FontSize',16);
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[1e-14 1]);
ylim(axes1,[1e-3 1e+6])
box(axes1,'on');
hold(axes1,'all');

for i=1:length(STD)
    
    loglog(B(:,5+STD(i))./N, p, COLOR{i}, 'LineWidth', SIZE(i), 'DisplayName',NAME{i});
    
end


% Create xlabel
xlabel('Mixing Ratio','FontSize',16);
% Create ylabel
ylabel('Pressure [Pa]','FontSize',16);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'EdgeColor',[1 1 1],'Location','BestOutside','FontSize',16);

set(gcf,'Position',[0,0,800,500]);