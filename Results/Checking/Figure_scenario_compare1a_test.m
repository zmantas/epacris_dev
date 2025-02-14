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

NAME = {'N_2','H_2O','CO','H','CH_4','NH_3','H_2S','H_2'};
    
STD = [55,7,20,3,21,9,45,53];

COLOR = {'k','b','r',[0.6 0.2 0],'c','m','g','y'};

% Create figure
figure1=figure;
axes1 = axes('Parent',figure1,'XScale','log','YScale','log','YDir','reverse','XMinorTick','on','YMinorTick','on',...
    'XTick',[1e-14 1e-12 1e-10 1e-8 1e-06 1e-4 1e-2 1],...
    'FontSize',16);
xlim(axes1,[1e-14 1]);
ylim(axes1,[1e-1 100e+5])
box(axes1,'on');
hold(axes1,'all');

for i=1:length(STD)
    
    loglog(B1(:,5+STD(i))./N1, p1, 'Color', COLOR{i}, 'LineWidth',3, 'DisplayName',NAME{i});
    
end

% Create label
xlabel('Mixing Ratio','FontSize',16);
ylabel('Pressure [Pa]','FontSize',16);
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest','FontSize',16);
set(gcf,'Position',[0,0,900,500]);

for i=1:length(STD)
    loglog(B2(:,5+STD(i))./N2, p2, 'Color', COLOR{i}, 'LineStyle', '--', 'LineWidth',3, 'DisplayName',NAME{i});
end

%
% for i=1:length(STD)
%     loglog(B3(:,5+STD(i))./N3, p, 'Color', COLOR{i}, 'LineStyle', ':', 'LineWidth',3, 'DisplayName',NAME{i});
% end

% saveas(gcf, 'Chem_CO2.fig');
% print -depsc2 'Chem_CO2.eps'
% close gcf