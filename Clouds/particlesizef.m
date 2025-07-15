% calculate particle size in exoplanet atmospheres

% input 
% g in SI
% T in K
% P in Pa
% M: mean molecular mass of the atmosphere; in g/mol
% MM: molecular mass of the condensable species; in g/mol
% KE: Eddy diffusion coefficient; in m2 s-1
% deltaP: difference between partial pressure and saturation vapor
% pressure, in Pa

% assume
% density of condensed material of water 1000 kg/m3
% accomodation factor of unity
% sig = 2

% output particle size in micron, and volumn in cm^3

function [r0,r1,r2,VP] = particlesizef(g,T,P,M,MM,KE,deltaP)

% Physical constant SI
HPLANCK = 6.626068E-34;
CLIGHT = 299792458;
KB = 1.3806503E-23;
AMU = 1.66053886E-27;
AVOGADRO = 6.022141E+23;
RGAS = 8.314472;
GRAVITY = 6.674E-11;

% Derived parameters
H = KB*T/M/AMU/g;
u = KE./H;
mu  = 8.76E-6*(293.85+72)./(293.85+72).*(T/293.85).^1.5; % SI
lambda = 2*mu./P./(8*M*1.0E-3/pi/8.314472./T).^0.5; % m
%KK = 4*KB*T/3./mu;
deltan=deltaP/KB/T;

% droplet
rho = 1.0E+3; % kg m-3
acc = 1;

% mass diffusion coefficient
D = 0.12E-4;

% Particle Size and Number    
Cc = 1;
fa = 1;
sig = 2;

% kappa=8*g*D*MM*AMU*deltan*exp(-2*log(sig)^2)/3/mu/u;
% disp(kappa);

for dump = 1:1e+3

cc = -(48*pi^2)^(1/3)*D*MM*AMU*fa*deltan/rho*exp(-log(sig)^2);
aa = rho*g/mu/(162*pi^2)^(1/3)/H*Cc*exp(-log(sig)^2);
bb = -u/H;
V = ((-bb+sqrt(bb^2-4*aa*cc))/2/aa)^(3/2);
d1 = (6*V/pi)^(1/3)*exp(-log(sig)^2);

kn  = lambda/d1;

Cc1  = 1+kn*(1.257+0.4*exp(-1.1/kn));
fa1  = (1+kn)/(1+2*kn*(1+kn)/acc);
if abs(Cc1-Cc)+abs(fa1-fa)<0.001
    Vs = V;
    break
else
    Cc = Cc1;
    fa = fa1;
end

end

r0 = (3*Vs/4/pi)^(1/3)*exp(-1.5*log(sig)^2)*1.0E+6;
r1 = (3*Vs/4/pi)^(1/3)*exp(-log(sig)^2)*1.0E+6;
r2 = (3*Vs/4/pi)^(1/3)*exp(-0.5*log(sig)^2)*1.0E+6;
VP = Vs*1.0E+6;

