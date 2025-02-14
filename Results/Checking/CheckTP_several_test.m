% Make Figure
clear all;
close all

c=jet(6);

% FILENAME = '../55cnce-iw4-test/NewTemperature_Jacob_converged1.dat';
FILENAME = '../55cnce/standard_200_test/IW+4/NewTemperature.dat';

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
    'XTick',[0 400 800 1200 1600 2000 2400 2800],...    
    'XMinorTick','on','YMinorTick','on', 'FontSize',16);
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 2800]);
ylim(axes1,[1e-3 2e+7])
box(axes1,'on');
hold(axes1,'all');
    
semilogy(t, p, 'Color', c(6,:), 'LineWidth', 2);

% Create xlabel
xlabel('Temperature [K]','FontSize',16);
% Create ylabel
ylabel('Pressure [Pa]','FontSize',16);
% Create legend
%legend1 = legend(axes1,'show');
%set(legend1,'EdgeColor',[1 1 1],'Location','BestOutside','FontSize',16);

set(gcf,'Position',[0,0,800,500]);



FILENAME = '../55cnce/standard_200_test/IW+2/NewTemperature.dat';

kb = 1.380658E-23;

% import data

IM=importdata(FILENAME,'\t');
z = IM(1:end,1);
t = IM(1:end,3);
p = 10.^IM(1:end,2);
N = p./t/kb*1.0E-6;

semilogy(t, p, 'Color', c(5,:), 'LineWidth', 2);

FILENAME = '../55cnce/standard_200_test/IW/NewTemperature.dat';

kb = 1.380658E-23;

% import data

IM=importdata(FILENAME,'\t');
z = IM(1:end,1);
t = IM(1:end,3);
p = 10.^IM(1:end,2);
N = p./t/kb*1.0E-6;

semilogy(t, p, 'Color', c(4,:), 'LineWidth', 2);

FILENAME = '../55cnce/standard_200_test/IW-2/NewTemperature.dat';

kb = 1.380658E-23;

% import data

IM=importdata(FILENAME,'\t');
z = IM(1:end,1);
t = IM(1:end,3);
p = 10.^IM(1:end,2);
N = p./t/kb*1.0E-6;

semilogy(t, p, 'Color', c(3,:), 'LineWidth', 2);

FILENAME = '../55cnce/standard_200_test/IW-4/NewTemperature.dat';

kb = 1.380658E-23;

% import data

IM=importdata(FILENAME,'\t');
z = IM(1:end,1);
t = IM(1:end,3);
p = 10.^IM(1:end,2);
N = p./t/kb*1.0E-6;

semilogy(t, p, 'Color', c(2,:), 'LineWidth', 2);

FILENAME = '../55cnce/standard_200_test/IW-6/NewTemperature.dat';

kb = 1.380658E-23;

% import data

IM=importdata(FILENAME,'\t');
z = IM(1:end,1);
t = IM(1:end,3);
p = 10.^IM(1:end,2);
N = p./t/kb*1.0E-6;

semilogy(t, p, 'Color', c(1,:), 'LineWidth', 2);













FILENAME = '../55cnce/standard_200_test_f03/IW+4/NewTemperature.dat';

kb = 1.380658E-23;

% import data

IM=importdata(FILENAME,'\t');
z = IM(1:end,1);
t = IM(1:end,3);
p = 10.^IM(1:end,2);
N = p./t/kb*1.0E-6;

semilogy(t, p, '--', 'Color', c(6,:), 'LineWidth', 2);


FILENAME = '../55cnce/standard_200_test_f03/IW+2/NewTemperature.dat';

kb = 1.380658E-23;

% import data

IM=importdata(FILENAME,'\t');
z = IM(1:end,1);
t = IM(1:end,3);
p = 10.^IM(1:end,2);
N = p./t/kb*1.0E-6;

semilogy(t, p, '--', 'Color', c(5,:), 'LineWidth', 2);

FILENAME = '../55cnce/standard_200_test_f03/IW/NewTemperature.dat';

kb = 1.380658E-23;

% import data

IM=importdata(FILENAME,'\t');
z = IM(1:end,1);
t = IM(1:end,3);
p = 10.^IM(1:end,2);
N = p./t/kb*1.0E-6;

semilogy(t, p, '--', 'Color', c(4,:), 'LineWidth', 2);

FILENAME = '../55cnce/standard_200_test_f03/IW-2/NewTemperature.dat';

kb = 1.380658E-23;

% import data

IM=importdata(FILENAME,'\t');
z = IM(1:end,1);
t = IM(1:end,3);
p = 10.^IM(1:end,2);
N = p./t/kb*1.0E-6;

semilogy(t, p, '--', 'Color', c(3,:), 'LineWidth', 2);

FILENAME = '../55cnce/standard_200_test_f03/IW-4/NewTemperature.dat';

kb = 1.380658E-23;

% import data

IM=importdata(FILENAME,'\t');
z = IM(1:end,1);
t = IM(1:end,3);
p = 10.^IM(1:end,2);
N = p./t/kb*1.0E-6;

semilogy(t, p, '--', 'Color', c(2,:), 'LineWidth', 2);

FILENAME = '../55cnce/standard_200_test_f03/IW-6/NewTemperature.dat';

kb = 1.380658E-23;

% import data

IM=importdata(FILENAME,'\t');
z = IM(1:end,1);
t = IM(1:end,3);
p = 10.^IM(1:end,2);
N = p./t/kb*1.0E-6;

semilogy(t, p, '--', 'Color', c(1,:), 'LineWidth', 2);