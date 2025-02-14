% Make Figure
%clear all;
%close all

% FILENAME = '../55cnce-iw4-test/NewTemperature_Jacob_converged1.dat';
FILENAME = '../k218b_soot_2/NewTemperature.dat';

kb = 1.380658E-23;

% import data

IM=importdata(FILENAME,'\t');
z = IM(1:end,1);
t = IM(1:end,3);
p = 10.^IM(1:end,2);
N = p./t/kb*1.0E-6;

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'XScale','lin','YScale','log','YDir','reverse',...
    'XTick',[0 100 200 300 400],...    
    'XMinorTick','on','YMinorTick','on', 'FontSize',16);
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 2800]);
ylim(axes1,[1e-3 1e+8])
box(axes1,'on');
hold(axes1,'all');
    
semilogy(t, p, 'k', 'LineWidth', 2);

% Create xlabel
xlabel('Temperature [K]','FontSize',16);
% Create ylabel
ylabel('Pressure [Pa]','FontSize',16);
% Create legend
%legend1 = legend(axes1,'show');
%set(legend1,'EdgeColor',[1 1 1],'Location','BestOutside','FontSize',16);

set(gcf,'Position',[0,0,800,500]);