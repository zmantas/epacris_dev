% Make Figure
clear all;
close all

FILENAME = '../55CncE_N2C/NewTemperature.dat';

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
    'XTick',[0 400 800 1200 1600 2000 2400 2800 3200 3600 4000 4400 4800],...    
    'XMinorTick','on','YMinorTick','on', 'FontSize',16);
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[1600 4000]);
ylim(axes1,[1e-0 1e+8])
box(axes1,'on');
hold(axes1,'all');
    
semilogy(t, p, 'b', 'LineWidth', 2,'DisplayName','N_2 Atmosphere');

% Create xlabel
xlabel('Temperature [K]','FontSize',16);
% Create ylabel
ylabel('Pressure [Pa]','FontSize',16);

set(gcf,'Position',[0,0,800,500]);

% import data
FILENAME = '../55CncE_O2/NewTemperature.dat';
IM=importdata(FILENAME,'\t');
z = IM(1:end,1);
t = IM(1:end,3);
p = 10.^IM(1:end,2);
N = p./t/kb*1.0E-6;

semilogy(t, p, 'r', 'LineWidth', 2,'DisplayName','O_2 Atmosphere');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'EdgeColor',[1 1 1],'Location','Best','FontSize',16);

t1=100:100:5000;
semilogy(t1, t1./t1*0.6E+5, '--k', 'LineWidth', 1);
semilogy(t1, t1./t1*3.4E+5, '--k', 'LineWidth', 1);

semilogy(t./t*2156, p, '--k', 'LineWidth', 1);
semilogy(t./t*2537, p, '--k', 'LineWidth', 1);