% Make Figure

global filename1;
global filename2;

KB=1.38064852E-23;

% import data

IM=importdata(filename1,'\t',2);
B1=IM.data;
z1 = B1(1:end,1);
t1 = B1(1:end,4);
p1 = B1(1:end,5);
N1 = p1./t1/KB*1e-6;

IM=importdata(filename2,'\t',2);
B2=IM.data;
z2 = B2(1:end,1);
t2 = B2(1:end,4);
p2 = B2(1:end,5);
N2 = p2./t2/KB*1e-6;

% select species

NAME = {'H','C','CH','CH_2','1CH_2','CH_3','C_2','C_2H','C_2H_3','C_2H_5'};

STD = [3,19,58,59,80,60,26,64,28,62];

COLOR = {[0.2 0.6 0],[0.6 0 0.6],'k','g','r','m','b','c',[0.6 0.2 0],[0 0.2 0.6]};

% Create figure
figure1=figure;
axes1 = axes('Parent',figure1,'XScale','log','YScale','log','YDir','reverse','XMinorTick','on','YMinorTick','on',...
    'XTick',[1e-14 1e-12 1e-10 1e-8 1e-06 1e-4 1e-2 1],...
    'FontSize',16);
xlim(axes1,[1e-20 1e-6]);
ylim(axes1,[1e-3 1e+4])
box(axes1,'on');
hold(axes1,'all');

for i=1:length(STD)
    
    loglog(B1(:,5+STD(i))./N1, p1/100, 'Color', COLOR{i}, 'LineWidth',3, 'DisplayName',NAME{i});
    
end

% Create label
xlabel('Mixing Ratio','FontSize',16);
ylabel('Pressure [mbar]','FontSize',16);
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest','FontSize',16);
set(gcf,'Position',[0,0,900,500]);

for i=1:length(STD)
    loglog(B2(:,5+STD(i))./N2, p2/100, 'Color', COLOR{i}, 'LineStyle', '--', 'LineWidth',3, 'DisplayName',NAME{i});
end

%
% for i=1:length(STD)
%     loglog(B3(:,5+STD(i))./N3, p, 'Color', COLOR{i}, 'LineStyle', ':', 'LineWidth',3, 'DisplayName',NAME{i});
% end

% saveas(gcf, 'Chem_CO2.fig');
% print -depsc2 'Chem_CO2.eps'
% close gcf