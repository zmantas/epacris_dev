% Make Figure
clear all;
close all

FILENAME1 = '../Earth_IW-6/NewTemperature.dat';
FILENAME1r = '../Earth_IW-6/NewTemperature_rc.dat';
FILENAME2 = '../Earth_IW-1/NewTemperature.dat';
FILENAME2r = '../Earth_IW-1/NewTemperature_rc.dat';
FILENAME3 = '../Earth_IW+4/NewTemperature.dat';
FILENAME3r = '../Earth_IW+4/NewTemperature_rc.dat';

kb = 1.380658E-23;

% import data

IM=importdata(FILENAME1,'\t');
z1 = IM(1:end,1);
t1 = IM(1:end,3);
p1 = 10.^IM(1:end,2);

IM=importdata(FILENAME1r,'\t');
z1r = IM(1:end,1);
t1r = IM(1:end,3);
p1r = 10.^IM(1:end,2);

IM=importdata(FILENAME2,'\t');
z2 = IM(1:end,1);
t2 = IM(1:end,3);
p2 = 10.^IM(1:end,2);

IM=importdata(FILENAME2r,'\t');
z2r = IM(1:end,1);
t2r = IM(1:end,3);
p2r = 10.^IM(1:end,2);

IM=importdata(FILENAME3,'\t');
z3 = IM(1:end,1);
t3 = IM(1:end,3);
p3 = 10.^IM(1:end,2);

IM=importdata(FILENAME3r,'\t');
z3r = IM(1:end,1);
t3r = IM(1:end,3);
p3r = 10.^IM(1:end,2);

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'XScale','lin','YScale','log','YDir','reverse',...
    'XTick',[0 400 800 1200 1600 2000 2400 2800],...    
    'XMinorTick','on','YMinorTick','on', 'FontSize',16);
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 2800]);
ylim(axes1,[1e-3 1e+8])
box(axes1,'on');
hold(axes1,'all');
    
semilogy(smooth(t1,9), p1, 'k', 'LineWidth', 2,'DisplayName','IW-6 Gray');
semilogy(smooth(t1r,9), p1r, '--k', 'LineWidth', 2,'DisplayName','IW-6 RC');
semilogy(smooth(t2,9), p2, 'b', 'LineWidth', 2,'DisplayName','IW-1 Gray');
semilogy(smooth(t2r,9), p2r, '--b', 'LineWidth', 2,'DisplayName','IW-1 RC');
semilogy(smooth(t3,9), p3, 'r', 'LineWidth', 2,'DisplayName','IW+4 Gray');
semilogy(smooth(t3r,9), p3r, '--r', 'LineWidth', 2,'DisplayName','IW+4 RC');

xlabel('Temperature [K]','FontSize',16);
ylabel('Pressure [Pa]','FontSize',16);
legend1 = legend(axes1,'show');
set(legend1,'EdgeColor',[1 1 1],'Location','Best','FontSize',16);

set(gcf,'Position',[0,0,800,500]);