% Make Figure

clear all;
close all

FILENAME = '../../Code/Result/Sun_0.1_SolarA_0/ConcentrationSTD_T.dat';

kb = 1.380658E-23;

% import data

IM=importdata(FILENAME,'\t',2);
B=IM.data;
N1 = sum(B(:,6:end)')';
z = B(1:end,1);
t = B(1:end,4);
p = B(1:end,5);
N = p./t/kb*1.0E-6;

He = (N-N1)./N;

mc = importdata('MoleCompo.dat');
mcH = mc(:,1);
mcO = mc(:,2);
mcN = mc(:,3);
mcC = mc(:,4);
mcS = mc(:,5);

HH = B(:,6:end)*mcH;
OO = B(:,6:end)*mcO;
NN = B(:,6:end)*mcN;
CC = B(:,6:end)*mcC;
SS = B(:,6:end)*mcS;

TT = HH + OO + NN + CC + SS;
var = [HH./TT, OO./TT, CC./TT, NN./TT, SS./TT];

% select species

NAME = {'H','O','C','N','S'};

SIZE = [3,3,3,3,3,3,3];
		 
COLOR = {'k','b','g','r','m'};

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'XScale','log','YScale','log','YDir','reverse',...
    'XTick',[1e-06 1e-4 1e-2 1],...    
    'XMinorTick','on','YMinorTick','on', 'FontSize',16);
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[1e-6 1]);
ylim(axes1,[1e-1 1e+8])
box(axes1,'on');
hold(axes1,'all');

for i=1:5
    
    loglog(var(:,i), p, COLOR{i}, 'LineWidth', SIZE(i), 'DisplayName',NAME{i});
    
end


% Create xlabel
xlabel('Mixing Ratio','FontSize',16);
% Create ylabel
ylabel('Pressure [Pa]','FontSize',16);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'EdgeColor',[1 1 1],'Location','BestOutside','FontSize',16);

set(gcf,'Position',[0,0,800,500]);