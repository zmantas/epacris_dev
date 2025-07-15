% cloud output
data = importdata('./CrossP/Cross_water_wavelength2.dat');
H2OL_r = data(:,1); % zero-order radius, in micron
H2OL_c = data(:,2:1388); % cross section per droplet, in cm2
data = importdata('./CrossP/Albedo_water_wavelength2.dat');
H2OL_a = data(:,2:1388);
data = importdata('./CrossP/Geo_water_wavelength2.dat');
H2OL_g = data(:,2:1388);

data = importdata('./CrossP/Cross_water_wavelength2.dat');
H2OI_r = data(:,1); % zero-order radius, in micron
H2OI_c = data(:,2:1388); % cross section per droplet, in cm2
data = importdata('./CrossP/Albedo_water_wavelength2.dat');
H2OI_a = data(:,2:1388);
data = importdata('./CrossP/Geo_water_wavelength2.dat');
H2OI_g = data(:,2:1388);

data = importdata('./CrossP/Cross_ammonia_wavelength2.dat');
NH3I_r = data(:,1);
NH3I_c = data(:,2:1388);
data = importdata('./CrossP/Albedo_ammonia_wavelength2.dat');
NH3I_a = data(:,2:1388);
data = importdata('./CrossP/Geo_ammonia_wavelength2.dat');
NH3I_g = data(:,2:1388);

% opacity
sig=2;
for j=1:length(zl)
    r2 = particlemsize(j);
    if cloudmden(j)<1e-12
        croa(j,1:1387)=0;
        alba(j,1:1387)=1;
        geoa(j,1:1387)=0;
    else
    r0 = r2*exp(-log(sig)^2); % micron
    VP = 4*pi/3*(r2*1.0E-6*exp(0.5*log(sig)^2))^3*1.0E+6*0.87; % g
    croa(j,:)= cloudmden(j)/VP*1.0E-3*10.^interp1(log10(NH3I_r),log10(NH3I_c),log10(max(0.01,min(r0,100)))); % cm-1
    alba(j,:)= interp1(log10(NH3I_r),NH3I_a,log10(max(0.01,min(r0,100))));
    geoa(j,:)= interp1(log10(NH3I_r),NH3I_g,log10(max(0.01,min(r0,100))));
    end
    r2 = particlesize(j);
    if cloudden(j)<1e-12
        crow(j,1:1387)=0;
        albw(j,1:1387)=1;
        geow(j,1:1387)=0;
    else
    r0 = r2*exp(-log(sig)^2);
    if tl(j)<273.16 % ice
        VP = 4*pi/3*(r2*1.0E-6*exp(0.5*log(sig)^2))^3*1.0E+6*0.92; % g
        crow(j,:)= cloudden(j)/VP*1.0E-3*10.^interp1(log10(H2OI_r),log10(H2OI_c),log10(max(0.01,min(r0,100)))); % cm-1
        albw(j,:)= interp1(log10(H2OI_r),H2OI_a,log10(max(0.01,min(r0,100))));
        geow(j,:)= interp1(log10(H2OI_r),H2OI_g,log10(max(0.01,min(r0,100))));
    else % liquid
        VP = 4*pi/3*(r2*1.0E-6*exp(0.5*log(sig)^2))^3*1.0E+6*1.0; % g
        crow(j,:)= cloudden(j)/VP*1.0E-3*10.^interp1(log10(H2OL_r),log10(H2OL_c),log10(max(0.01,min(r0,100)))); % cm-1
        albw(j,:)= interp1(log10(H2OL_r),H2OL_a,log10(max(0.01,min(r0,100))));
        geow(j,:)= interp1(log10(H2OL_r),H2OL_g,log10(max(0.01,min(r0,100))));
    end
    end
end

% dlmwrite([outdir,'Cross_H2O_R.dat'],crow,'delimiter','\t','precision','%.6e');
% dlmwrite([outdir,'albedo_H2O_R.dat'],albw,'delimiter','\t','precision','%.6e');
% dlmwrite([outdir,'geo_H2O_R.dat'],geow,'delimiter','\t','precision','%.6e');
% dlmwrite([outdir,'Cross_NH3_R.dat'],croa,'delimiter','\t','precision','%.6e');
% dlmwrite([outdir,'albedo_NH3_R.dat'],alba,'delimiter','\t','precision','%.6e');
% dlmwrite([outdir,'geo_NH3_R.dat'],geoa,'delimiter','\t','precision','%.6e');

% dlmwrite([outdir,'Cross_H2O_C.dat'],crow(idbottom:end,:),'delimiter','\t','precision','%.6e');
% dlmwrite([outdir,'albedo_H2O_C.dat'],albw(idbottom:end,:),'delimiter','\t','precision','%.6e');
% dlmwrite([outdir,'geo_H2O_C.dat'],geow(idbottom:end,:),'delimiter','\t','precision','%.6e');
% dlmwrite([outdir,'Cross_NH3_C.dat'],croa(idbottom:end,:),'delimiter','\t','precision','%.6e');
% dlmwrite([outdir,'albedo_NH3_C.dat'],alba(idbottom:end,:),'delimiter','\t','precision','%.6e');
% dlmwrite([outdir,'geo_NH3_C.dat'],geoa(idbottom:end,:),'delimiter','\t','precision','%.6e');

% full cloud information
f = fopen([outdir,'Cloudinfo.dat'],'w');
fprintf(f,'%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t','z','z0','z1','T','P');
fprintf(f,'%s\t\t','H2O Den');
fprintf(f,'%s\t\t','NH3 Den');
fprintf(f,'%s\t\t','H2O Rad');
fprintf(f,'%s\t\t','NH3 Rad');
fprintf(f,'\n');
fprintf(f,'%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t','km','km','km','K','Pa');
fprintf(f,'%s\t\t','g/L');
fprintf(f,'%s\t\t','g/L');
fprintf(f,'%s\t\t','micron');
fprintf(f,'%s\t\t','micron'); % radius, 2-moment (surface average) for lognormal with sigma=2
fprintf(f,'\n');
for j=1:length(zl)
    fprintf(f,'%.6f\t%.6f\t%.6f\t%.6f\t%.6e\t',zl(j),z0(j),z1(j),tl(j),pl(j));
    fprintf(f,'%.6e\t',cloudden(j));
    fprintf(f,'%.6e\t',cloudmden(j));
    fprintf(f,'%.6e\t',particlesize(j));
    fprintf(f,'%.6e\t',particlemsize(j));
    fprintf(f,'\n');
end
fclose(f);

