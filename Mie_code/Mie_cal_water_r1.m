% calculate the Mie theory scattering parameters to be included in the Main
% C loop

% Hypothetical AER with n=1.5 (NaCl, KCl), and k=0

% Using the Lognormal Distribution of particle sizes

clear all
close all

% SET THE WAVELENGTH OF THE CALCULATION
%lam1 = 10.^(-1:0.01:2); % micron
NLAMBDA = 16000;
LAMBDALOW = 0.1*1.e-6;
LAMBDAHIGH = 200*1.e-6;
lam1= zeros(NLAMBDA, 1);
start = log10(LAMBDALOW);
interval = log10(LAMBDAHIGH) - log10(LAMBDALOW);
dlambda = interval / (NLAMBDA - 1);

for i = 0:(NLAMBDA-1)
    lam1(i+1) = 10.0^(start + i * dlambda) * 1.e+6;
end

% LOAD THE REFRACTIVE INDICES
dum=importdata('H2OLiquid_N','\t',2);
refidx=dum.data;

% SET THE PARTICLE DIAMETER OF THE CALCULATION
d = 10.^(-3:0.25:2); % micron

for ii=1:length(d)

% SET THE RADIUS OF THE PARTICLE
a = d(ii)/2;
% SET THE STANDARD DEVIATION OF PARTICLE SIZE
sig = 2; % lognormal dispersion
a0 = a*exp(-log(sig)^2);
r = a0*sig.^(-5:0.01:5);
lnr = log(r);
fr = 1/log(sig)/sqrt(2*pi)*exp(-(lnr-log(a0)).^2/2/log(sig)^2);

% Calculate from the Mie theory

for i=1:length(lam1)
    
    lambda = lam1(i);
    nl(i) = interp1(refidx(:,1),refidx(:,2),lambda);
    kl(i) = interp1(refidx(:,1),refidx(:,3),lambda);
    m(i)=nl(i)+kl(i)*sqrt(-1);
    
    clear fext fsca fg
    for j=1:length(r)
        x=2*pi*r(j)/lambda;
        result = Mie(m(i),x);
        fext(j) = result(4)*pi*r(j)^2*1E-8; % in cm2
        fsca(j) = result(5)*pi*r(j)^2*1E-8; % in cm2
        fg(j)   = result(8)*fsca(j);
        if isnan(fg(j))
            fg(j) = 0;
        end
    end
    
    sigext(i,ii) = trapz(lnr,fext.*fr); % in cm2
    sigsca(i,ii) = trapz(lnr,fsca.*fr); % in cm2
    g(i,ii) = trapz(lnr,fg.*fr)/sigsca(i,ii);
    
end

disp(sprintf('%d/%d, d=%d',ii, length(d),d(ii)));

end

foutc = fopen('result/Cross_water_110123.dat','w');
fouta = fopen('result/Albedo_water_110123.dat','w');
foutg = fopen('result/Geo_water_110123.dat','w');

for i=1:length(lam1)
    fprintf(foutc, '%.6e\t', lam1(i));
    for j=1:length(d)
        fprintf(foutc, '%.6e\t', sigext(i,j));
    end
    fprintf(foutc, '\n');
    fprintf(fouta, '%.6e\t', lam1(i));
    for j=1:length(d)
        fprintf(fouta, '%.6e\t', sigsca(i,j)/sigext(i,j));
    end
    fprintf(fouta, '\n');
    fprintf(foutg, '%.6e\t', lam1(i));
    for j=1:length(d)
        fprintf(foutg, '%.6f\t', g(i,j));
    end
    fprintf(foutg, '\n');
end

fclose(foutc);
fclose(fouta);
fclose(foutg);
