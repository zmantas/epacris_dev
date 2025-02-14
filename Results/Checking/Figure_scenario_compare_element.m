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
N1r = sum(B1(:,6:end),2);

IM=importdata(filename2,'\t',2);
B2=IM.data;
z2 = B2(1:end,1);
t2 = B2(1:end,4);
p2 = B2(1:end,5);
N2 = p2./t2/KB*1e-6;
N2r = sum(B2(:,6:end),2);

% elemental composition
compo = importdata('SpeciesName.txt');
compo = compo.data;

compo_H = compo(:,1);
H1 = B1(:,6:end)*compo_H;
H2 = B2(:,6:end)*compo_H;

compo_O = compo(:,2);
O1 = B1(:,6:end)*compo_O;
O2 = B2(:,6:end)*compo_O;

compo_N = compo(:,3);
N1 = B1(:,6:end)*compo_N;
N2 = B2(:,6:end)*compo_N;

compo_C = compo(:,4);
C1 = B1(:,6:end)*compo_C;
C2 = B2(:,6:end)*compo_C;

compo_S = compo(:,5);
S1 = B1(:,6:end)*compo_S;
S2 = B2(:,6:end)*compo_S;

COLOR = {'k','b','r','c','m'};
NAME = {'H','O','N','C','S'};

% Create figure
figure1=figure;
axes1 = axes('Parent',figure1,'XScale','log','YScale','log','YDir','reverse','XMinorTick','on','YMinorTick','on',...
    'XTick',[1e-14 1e-12 1e-10 1e-8 1e-06 1e-4 1e-2 1],...
    'FontSize',16);
%xlim(axes1,[1e-6 1]);
ylim(axes1,[1e-1 100e+5])
box(axes1,'on');
hold(axes1,'all');
   
loglog(H1./N1r, p1, 'Color', COLOR{1}, 'LineWidth',3, 'DisplayName',NAME{1});
loglog(O1./N1r, p1, 'Color', COLOR{2}, 'LineWidth',3, 'DisplayName',NAME{2});
loglog(N1./N1r, p1, 'Color', COLOR{3}, 'LineWidth',3, 'DisplayName',NAME{3});
loglog(C1./N1r, p1, 'Color', COLOR{4}, 'LineWidth',3, 'DisplayName',NAME{4});
loglog(S1./N1r, p1, 'Color', COLOR{5}, 'LineWidth',3, 'DisplayName',NAME{5});

% Create label
xlabel('Mixing Ratio','FontSize',16);
ylabel('Pressure [Pa]','FontSize',16);
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest','FontSize',16);
set(gcf,'Position',[0,0,900,500]);

loglog(H2./N2r, p2, 'Color', COLOR{1},'LineStyle', '--', 'LineWidth',3, 'DisplayName',NAME{1});
loglog(O2./N2r, p2, 'Color', COLOR{2},'LineStyle', '--', 'LineWidth',3, 'DisplayName',NAME{2});
loglog(N2./N2r, p2, 'Color', COLOR{3},'LineStyle', '--', 'LineWidth',3, 'DisplayName',NAME{3});
loglog(C2./N2r, p2, 'Color', COLOR{4},'LineStyle', '--', 'LineWidth',3, 'DisplayName',NAME{4});
loglog(S2./N2r, p2, 'Color', COLOR{5},'LineStyle', '--', 'LineWidth',3, 'DisplayName',NAME{5});
%
% for i=1:length(STD)
%     loglog(B3(:,5+STD(i))./N3, p, 'Color', COLOR{i}, 'LineStyle', ':', 'LineWidth',3, 'DisplayName',NAME{i});
% end

% saveas(gcf, 'Chem_CO2.fig');
% print -depsc2 'Chem_CO2.eps'
% close gcf