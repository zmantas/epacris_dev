% Temperature profile of grey atmosphere
% Output T in K and P in Pascal
% LoopVar: variation of successive loop to quality check
% watermix: mixing ratio of water vapor
% cloudden: density of cloud, in g/L
% particlesize: quadratic mean of particle size, in micron
% cloudopacity: opacity of cloud particles in cm-1
% cloudcompo: NH3/H2O molar ratio in water cloud
% albedoM: Bond albedo in 5 optical bands
% albedo: Bond albedo
% Input
% Surface Gravity -- g in SI unit
% Metallicity -- meta [M/H] -- assuming linear dependency of metallicity
    % and opacities
% Irradiation Temperature -- Tirr in K
% Intrinsic Temperature -- Tint in K
% Eddy Diffusion Coefficient -- KE in m2/s
% Formulation by Guillot (2010)
% Assuming full heat redistribution
% Gamma value recommended by Guillot (2010), scaled so at 600 K give 0.07
% Check atmosphere stability of H2-He Mix
% Check Water Condensation
% Use dry adiabatic lapse rate when condensation dost not occur
% Use moist adiabatic lapse rate when condenation occur
% Particle cross section calculated by Mie theory assuming sigma=2 lognrmal
% assuming 1 g/cm3 for water
% Produce standard output for RT calculation
% The opar is to correct the Rossland mean opacity at low pressure
% Seperate output from 10x tropopause pressure in 70 layers for
% photochemical calculations
% Include NH4SH cloud formation

function [ T, P, LoopVar, watermix, cloudden, particlesize, cloudopacity, cloudcompo, ammoniamix, cloudmden, particlemsize, cloudmopacity, albedoD, albedo, Pref] ...
    = tpgrey_albedoD_gasplanetH5r_photo( g, meta, opar, xH2, Tirr, Tint, KE, cloudfrac, methanemix, hazemix, DIRName )

    global ha ha_m ha_t;

    % Physical Constant
    kb = 1.3806488E-23; % Boltzmann constant in SI
    AMU = 1.66053892E-27; % Atomic mass unit in SI
    NA=6.0221413E+23; % Avogadro's number
    SB=5.670374419e-8; % Stefan-Bolzmann constant
    ATM=101325;
    
    % Wavelength
    wave = 100:1:5000; % nm

    % Solar Spectrum
    data = importdata('RayleightOpa/solar.txt');
    solar = interp1(data(:,1),data(:,2),wave)';

    % Methane Opacity
    data = importdata('MethaneOpa/CH4cross.txt');
    crossCH4 = interp1(data(:,1),data(:,2),wave);

    % Ammonia Opacity
    data = importdata('MethaneOpa/NH3cross.txt');
    crossNH3 = interp1(data(:,1),data(:,2),wave);

    % Water Opacity
    data = importdata('MethaneOpa/H2Ocross.txt');
    crossH2O = interp1(data(:,1),data(:,2),wave);
    
    % Rayleigh Opacity
    x = wave/1000.0;
    nHe = 1 + 0.01470091./(423.98-x.^(-2));
    nH2 = 1.0001393*x./x;
    DenS=101325.0/1.3806488E-23/273.0*1.0E-6;
    crossHe=8.0*pi^3*(nHe.^2-1).^2/3./(wave*1.0E-7).^4/DenS/DenS;
    crossH2=8.0*pi^3*(nH2.^2-1).^2/3./(wave*1.0E-7).^4/DenS/DenS;
    crossRay= crossH2*xH2 + crossHe*(1-xH2);

    % Set up pressure grid
    P = 10.^(-4:0.02:9); % Pascal
    
    % Convert pressure to column mass g/cm2
    m = P/g*1.0E-1;
    
    % Initial Temperature
    T = Tirr*P./P;
    
    % Atmosphere Composition
    [ fH2, fHe, fH2O, fCH4, fNH3, fH2S, MMM ] = metacompo( xH2, meta );
    % fCH4=methanemix;
    % fH2O=fH2O*2;
    % fNH3=fNH3*2;

    % disp(fCH4);
    
    % Load Mie Calculation Results
    data = importdata('CrossP/cross_H2OLiquid_M.dat');
    H2OL_r = data(:,1); % zero-order radius, in micron
    H2OL_c = data(:,4); % cross section per droplet, in cm2
    data = importdata('CrossP/cross_H2OIce_M.dat');
    H2OI_r = data(:,1);  
    H2OI_c = data(:,4);
    data = importdata('CrossP/cross_NH3Ice_M.dat');
    NH3I_r = data(:,1);  
    NH3I_c = data(:,4);
    
    % miu value integral
    %miu = [0.0001,0.1429,0.2857,0.4286,0.5714,0.7143,0.8571,1.0000]';
    miu = linspace(0.0001,0.9999,10)';
    miuc= trapz(miu,miu.^2)/(1/3);
    imiu= 1./miu;
    
    % Loop to update temperature
    LoopMax = 10;
    LoopCri = 1E-3;
    LoopVar = 1;
    LoopID  = 1;
    E2 = importdata('E2.dat');
    albedo = 0.343; % initial albedo
    
%     particlemsize = 1e-16;
%     particlesize  = 1e-16;
%     sig = 2;
 
    while LoopID<LoopMax && LoopVar>LoopCri
        
        % Calculate Optical Depth
        data=importdata('MeanOpacity/MH0-RM.txt');
        tg=data(1,2:end);
        pg=data(2:end,1);
        rs=data(2:end,2:end)*opar;
        for i=1:length(P)
            if i==1
                opa = 10.^interp2(tg,log10(pg),log10(rs),min(max(75,T(i)),4000),log10(min(3E+8,max(3E+2,P(i)))));
                tau(i)=m(i)*opa*5^meta;
            else
                opa = 10.^interp2(tg,log10(pg),log10(rs),min(max(75,T(i)),4000),log10(min(3E+8,max(3E+2,P(i)))));
                tau(i)=tau(i-1)+(m(i)-m(i-1))*opa*5^meta;
            end
        end
        
        % Calculate Gamma
        gamma = 0.13*sqrt(Tirr/2000);
        
        % Calculate Tnew
        Teq=Tirr*(1/2)^(1/2);
        EE = interp1(E2(:,1),E2(:,2),min(1E+5,max(1E-3,gamma*tau)));
        Tnew4 = 3*Tint^4/4*(2/3+tau) + (1-albedo)*3*Teq^4/4*(2/3+2/3/gamma*(1+(gamma*tau/2-1).*exp(-gamma*tau))+2*gamma/3*(1-tau.^2/2).*EE);
        Tnew = Tnew4.^0.25;
        
        % Calculate Molar Heat Capacity
        for i=1:length(Tnew)
            th=min(6000,max(298,Tnew(i)));
            t=th/1000.0;
            % H2
            if th<1000.0
                A=33.066178;
                B=-11.363417;
                C=11.432816;
                D=-2.772874;
                E=-0.158558;
            elseif th<2500.0
                A=18.563083;
                B=12.257357;
                C=-2.859786;
                D=0.268238;
                E=1.97799;
            else
                A=43.41356;
                B=-4.293079;
                C=1.272428;
                D=-0.096876;
                E=-20.533862;
            end
            cpH2(i) = A + B*t +  C*t*t + D*t*t*t + E/t/t;
            % He
            	A=20.78603;
                B=4.850638E-10;
                C=-1.582916E-10;
                D=1.525102E-11;
                E=3.196347E-11;
            cpHe(i) = A + B*t +  C*t*t + D*t*t*t + E/t/t;
            % CH4
            if th<1300.0
                A=-0.703029;
                B=108.4773;
                C=-42.52157;
                D=5.862788;
                E=0.678565;
            else
                A=85.81217;
                B=11.26467;
                C=-2.114146;
                D=0.138190;
                E=-26.42221;
            end
            cpCH4(i) = A + B*t +  C*t*t + D*t*t*t + E/t/t;
            % H2O
            th=min(6000,max(500,Tnew(i)));
            t=th/1000.0;
            if th<1700.0
                A=30.09200;
                B=6.832514;
                C=6.793435;
                D=-2.534480;
                E=0.082139;
            else
                A=41.96426;
                B=8.622053;
                C=-1.499780;
                D=0.098119;
                E=-11.15764;
            end
            cpH2O(i) = A + B*t +  C*t*t + D*t*t*t + E/t/t;
            
        end
        cp = cpH2*fH2+cpHe*fHe+cpH2O*fH2O+cpCH4*fCH4;
        cp = (cpH2*fH2+cpHe*fHe)/(fH2+fHe);
        
        % Calculate lapse rate
        cpm = cp/(MMM*0.001); % convert to J/kg/K;
        lapse = kb/MMM/AMU./cpm; % dry adiabatic lapse rate for dlnT/dlnP dimensionless
        
        % Process the temperature profile
        for i=1:length(P)-1
            dlnTdlnP=(log(Tnew(i+1))-log(Tnew(i)))/(log(P(i+1))-log(P(i)));
            if dlnTdlnP>lapse(i)
               Tnew(i+1)=exp(log(Tnew(i))+(log(P(i+1))-log(P(i)))*lapse(i));
            end
        end
        
        % Check atmosphere stability and cloud formation
        % WaterSolar = 8.8E-4; % Water mixing ratio at Solar abundances
        Lw1 = 2500800; % latent heat of water evaporation, J/kg at 0C
        Lw2 = 2834100; % latent heat of water sublimation, J/kg at 0C
        LNH3= 1371200; % latent heat of NH3 vaporization, J/kg at -33.5C
        RGAS = 8.3144621; % Universal Gas Constant, SI
        watermix = fH2O*P./P;
        ammoniamix = fNH3*P./P;
        sulfidemix = fH2S*P./P;
        for i=length(P):-1:2
            ifwaterc = Tnew(i)<647.096 && waterpressure(Tnew(i))<=P(i)*watermix(i);
            ifammoniac = Tnew(i)<300.0 && ammoniapressure(Tnew(i))<=P(i)*ammoniamix(i);
            
            if ifwaterc && ifammoniac==0 % condensation of water only
                deltaP = P(i)*watermix(i)-waterpressure(Tnew(i));
                watermixn= waterpressure(Tnew(i))/P(i);
                cloudden(i)= max((watermix(i)-watermixn)*0.018.*P(i)./RGAS./Tnew(i),1e-36); % kg/m^3, g/L
                watermix(i)= watermixn;
                if Tnew(i)>=273.16
                    [ammoniamix(i),cloudcompo(i)]=ammoniahenry_var(Tnew(i),P(i),cloudden(i),ammoniamix(i));
                end
                cloudden(i)=cloudden(i)+watermixn*cloudcompo(i)*0.018.*P(i)./RGAS./Tnew(i); % update cloudden due to ammonia dissolution
                % disp([cloudcompo(i)]);
                % calculate cloud particle size
                [r0,r1,r2,VP] = particlesizef(g,Tnew(i),P(i),MMM,18.0,KE,deltaP);
                particlesize(i)=r2;
                % calculate moist lapse rate and cloud opacity profile
                if Tnew(i)<273.16 % ice
                    lapsem=lapse(i)*(1+Lw2*0.018*watermix(i)/RGAS/Tnew(i))/(1+Lw2^2*0.018^2*watermix(i)/MMM/0.001/RGAS/Tnew(i)^2/cpm(i));
                    % disp([lapse(i),lapsem]);
                    cloudopacity(i)= cloudden(i)/VP*1.0E-3*10.^interp1(log10(H2OI_r),log10(H2OI_c),log10(max(0.01,min(r0,100)))); % cm-1
                else % liquid
                    lapsem=lapse(i)*(1+Lw1*0.018*watermix(i)/RGAS/Tnew(i))/(1+Lw1^2*0.018^2*watermix(i)/MMM/0.001/RGAS/Tnew(i)^2/cpm(i));
                    % disp([lapse(i),lapsem]);
                    cloudopacity(i)= cloudden(i)/VP*1.0E-3*10.^interp1(log10(H2OL_r),log10(H2OL_c),log10(max(0.01,min(r0,100)))); % cm-1
                end
                % see if convect
                dlnTdlnP=(log(Tnew(i-1))-log(Tnew(i)))/(log(P(i-1))-log(P(i)));
                if dlnTdlnP>lapsem
                    %disp([log10(P(i)),lapse(i),dlnTdlnP]);
                    Tnew(i-1)=exp(log(Tnew(i))+(log(P(i-1))-log(P(i)))*lapsem);
                end
                % water vapor decrease
                % watermixn = waterpressure(Tnew(i-1))/P(i-1);
                watermix(i-1) = watermix(i);
                ammoniamix(i-1) = ammoniamix(i);
                cloudmden(i) = 1e-36;
                particlemsize(i) = 1e-36;
                cloudmopacity(i) = 1e-36;
            end
            
            if ifammoniac && ifwaterc==0 % condensation of ammonia only 
                deltaP = ammoniamix(i)*P(i)-ammoniapressure(Tnew(i));
                ammoniamixn= ammoniapressure(Tnew(i))/P(i);
                cloudmden(i)= max((ammoniamix(i)-ammoniamixn)*0.017.*P(i)./RGAS./Tnew(i),1e-36); % kg/m^3, g/L
                ammoniamix(i)= ammoniamixn;
                % calculate cloud particle size
                [r0,r1,r2,VP] = particlesizef(g,Tnew(i),P(i),MMM,17.0,KE,deltaP);
                particlemsize(i)=r2;
                % calculate moist lapse rate and cloud opacity profile
                lapsem=lapse(i)*(1+LNH3*0.017*ammoniamix(i)/RGAS/Tnew(i))/(1+LNH3^2*0.017^2*ammoniamix(i)/MMM/0.001/RGAS/Tnew(i)^2/cpm(i));
                cloudmopacity(i)= cloudmden(i)/VP*1.0E-3*10.^interp1(log10(NH3I_r),log10(NH3I_c),log10(max(0.01,min(r0,100)))); % cm-1
                % see if convect
                dlnTdlnP=(log(Tnew(i-1))-log(Tnew(i)))/(log(P(i-1))-log(P(i)));
                if dlnTdlnP>lapsem
                    Tnew(i-1)=exp(log(Tnew(i))+(log(P(i-1))-log(P(i)))*lapsem);
                end
                watermix(i-1) = watermix(i);
                ammoniamix(i-1) = ammoniamix(i);
                cloudden(i) = 1e-36;
                particlesize(i) = 1e-36;
                cloudopacity(i) = 1e-36;
                cloudcompo(i) = 0;
            end
            
            if ifammoniac && ifwaterc % both water and ammonia condense
                deltaP1 = P(i)*watermix(i)-waterpressure(Tnew(i));
                watermixn= waterpressure(Tnew(i))/P(i);
                cloudden(i)= max((watermix(i)-watermixn)*0.018.*P(i)./RGAS./Tnew(i),1e-36); % kg/m^3, g/L
                watermix(i)= watermixn;
                deltaP2 = P(i)*ammoniamix(i)-ammoniapressure(Tnew(i));
                ammoniamixn= ammoniapressure(Tnew(i))/P(i);
                cloudmden(i)= max((ammoniamix(i)-ammoniamixn)*0.017.*P(i)./RGAS./Tnew(i),1e-36); % kg/m^3, g/L
                ammoniamix(i)= ammoniamixn;
                % calculate cloud particle size and cloud opacity profile
                [r0,r1,r2,VP] = particlesizef(g,Tnew(i),P(i),MMM,18.0,KE,deltaP1);
                particlesize(i) = r2;
                if Tnew(i)<273.16 % ice
                    % lapsem=lapse(i)*(1+Lw2*0.018*watermix(i)/RGAS/Tnew(i))/(1+Lw2^2*0.018^2*watermix(i)/MMM/0.001/RGAS/Tnew(i)^2/cpm(i));
                    % disp([lapse(i),lapsem]);
                    cloudopacity(i)= cloudden(i)/VP*1.0E-3*10.^interp1(log10(H2OI_r),log10(H2OI_c),log10(max(0.01,min(r0,100)))); % cm-1
                else % liquid
                    % lapsem=lapse(i)*(1+Lw1*0.018*watermix(i)/RGAS/Tnew(i))/(1+Lw1^2*0.018^2*watermix(i)/MMM/0.001/RGAS/Tnew(i)^2/cpm(i));
                    % disp([lapse(i),lapsem]);
                    cloudopacity(i)= cloudden(i)/VP*1.0E-3*10.^interp1(log10(H2OL_r),log10(H2OL_c),log10(max(0.01,min(r0,100)))); % cm-1
                end
                % calculate cloud particle size
                [r0,r1,r2,VP] = particlesizef(g,Tnew(i),P(i),MMM,17.0,KE,deltaP2);
                particlemsize(i) = r2;
                cloudmopacity(i)= cloudmden(i)/VP*1.0E-3*10.^interp1(log10(NH3I_r),log10(NH3I_c),log10(max(0.01,min(r0,100)))); % cm-1
                % calculate moist lapse rate
                if Tnew(i)<273.16 % ice
                    lapsem=lapse(i)*(1+Lw2*0.018*watermix(i)/RGAS/Tnew(i)+LNH3*0.017*ammoniamix(i)/RGAS/Tnew(i))/(1+Lw2^2*0.018^2*watermix(i)/MMM/0.001/RGAS/Tnew(i)^2/cpm(i)+LNH3^2*0.017^2*ammoniamix(i)/MMM/0.001/RGAS/Tnew(i)^2/cpm(i));
                else % liquid
                    lapsem=lapse(i)*(1+Lw1*0.018*watermix(i)/RGAS/Tnew(i)+LNH3*0.017*ammoniamix(i)/RGAS/Tnew(i))/(1+Lw1^2*0.018^2*watermix(i)/MMM/0.001/RGAS/Tnew(i)^2/cpm(i)+LNH3^2*0.017^2*ammoniamix(i)/MMM/0.001/RGAS/Tnew(i)^2/cpm(i));
                end
                %disp([P(i),lapsem/lapse(i)]);
                % see if convect
                dlnTdlnP=(log(Tnew(i-1))-log(Tnew(i)))/(log(P(i-1))-log(P(i)));
                if dlnTdlnP>lapsem
                    Tnew(i-1)=exp(log(Tnew(i))+(log(P(i-1))-log(P(i)))*lapsem);
                end
                watermix(i-1) = watermix(i);
                ammoniamix(i-1) = ammoniamix(i);
                cloudcompo(i) = 0;
            end
            
            if  ifammoniac==0 && ifwaterc==0% no condensation
                cloudden(i) = 1e-36;
                particlesize(i) = 1e-36;
                cloudopacity(i) = 1e-36;
                cloudcompo(i) = 0;
                cloudmden(i) = 1e-36;
                particlemsize(i) = 1e-36;
                cloudmopacity(i) = 1e-36;
                dlnTdlnP=(log(Tnew(i-1))-log(Tnew(i)))/(log(P(i-1))-log(P(i)));
                if dlnTdlnP>lapse(i)
                    %disp([log10(P(i)),lapse(i),dlnTdlnP]);
                    Tnew(i-1)=exp(log(Tnew(i))+(log(P(i-1))-log(P(i)))*lapse(i));
                end
                watermix(i-1) = watermix(i);
                ammoniamix(i-1) = ammoniamix(i);
            end
            
            % NH4SH
            pammonia=P(i)*ammoniamix(i)/ATM;
            psulfide=P(i)*sulfidemix(i)/ATM;
            eqnh4sh =10^(14.82-4705/Tnew(i));
            ifnh4sh = (pammonia*psulfide)>eqnh4sh; % compare with equilibrium constant
            if ifnh4sh
                pdelta = (pammonia+psulfide-sqrt((pammonia+psulfide)^2-4*(pammonia*psulfide-eqnh4sh)))/2;
                ammoniamixn = (pammonia-pdelta)*ATM/P(i);
                sulfidemixn = (psulfide-pdelta)*ATM/P(i);
                % count as NH3 clouds
                cloudmden(i)= cloudmden(i) + max((ammoniamix(i)-ammoniamixn)*0.051.*P(i)./RGAS./Tnew(i),1e-36); % kg/m^3, g/L
                ammoniamix(i)= ammoniamixn;
                sulfidemix(i)= sulfidemixn;
                % calculate cloud particle size if not calculated
                if particlemsize(i) < 1e-20
                    [r0,r1,r2,VP] = particlesizef(g,Tnew(i),P(i),MMM,51.0,KE,pdelta*ATM);
                    particlemsize(i)=r2;
                end
                watermix(i-1) = watermix(i);
                ammoniamix(i-1) = ammoniamix(i);
                sulfidemix(i-1) = sulfidemix(i);
            else
                sulfidemix(i-1) = sulfidemix(i);
            end
            
        end
        
        % Update albedo
        i=1;
        refdepth=0;
        wrefdepth=0;
        clouddepth=0;
        wclouddepth=0;
        clouddepth1=1e-36;
        wclouddepth1=1e-36;
        rayleighdepth=0;
        albedod = (miu./miu)*(wave./wave);
        albedod0= (miu./miu)*(wave./wave);
        cloudmass = 0;
        cloudmmass = 0;
        while i<length(P)
            sh=kb*T(i)/MMM/AMU/g*100; % cm
            clouddepth=clouddepth+(cloudopacity(i)+cloudmopacity(i))*log(P(i+1)/P(i))*sh; %cloud
            wclouddepth=wclouddepth+cloudopacity(i)*log(P(i+1)/P(i))*sh; % water cloud
            clouddepth1(i+1)=clouddepth1(i)+(cloudopacity(i+1)+cloudmopacity(i+1))*log(P(i+1)/P(i))*sh; %cloud
            wclouddepth1(i+1)=wclouddepth1(i)+cloudopacity(i+1)*log(P(i+1)/P(i))*sh; % water cloud
            cloudmass = cloudmass+cloudden(i)*log(P(i+1)/P(i))*sh/100;
            cloudmmass = cloudmmass+cloudmden(i)*log(P(i+1)/P(i))*sh/100;
            rayleighdepth=rayleighdepth+crossRay*(P(i+1)-P(i))/g*1.0E-4/MMM/AMU;
            refdepth=clouddepth+rayleighdepth;
            wrefdepth=wclouddepth+rayleighdepth;
            ridx=(imiu*refdepth<1);
            wridx=(imiu*wrefdepth<1);
            gasopa=(crossCH4*fCH4+crossH2O*watermix(i)+crossNH3*ammoniamix(i))*(P(i+1)-P(i))/g*1.0E-4/MMM/AMU;
            albedod=albedod.*exp(-2*ridx.*(imiu*gasopa));
            albedod0=albedod0.*exp(-2*wridx.*(imiu*gasopa));
            redx(i,:)=1-ridx(10,:);
            wedx(i,:)=1-wridx(10,:);
            i=i+1;
        end
        
        % semilogx(wave,albedod0);
        
        albedoD = 2*trapz(miu,albedod.*(miu*(wave./wave)).^2)/miuc;
        albedoD0 = 2*trapz(miu,albedod0.*(miu*(wave./wave)).^2)/miuc;

        %[id1,id2]=max(redx);
        %Pref=P(id2);
        %[id1,id2]=max(wedx);
        %wPref=P(id2);
        %plot(wave,Pref);
        %Pref=1;
        [clouddepth1,cuidx]=unique(clouddepth1);
        Pref=10^interp1(log10(clouddepth1),log10(P(cuidx)),0);
        [wclouddepth1,wcuidx]=unique(wclouddepth1);
        wPref=10^interp1(log10(wclouddepth1),log10(P(wcuidx)),0);

        disp([cloudmass,cloudmmass]); % kg m-2

        disp([clouddepth,wclouddepth]);

%         % Update albedo
%         i=1;
%         clouddepth=0;
%         gasdepth=0;
%         rayleighdepth=0;
%         while min(clouddepth)<1 && i<length(P)
%             sh=kb*T(i)/MMM/AMU/g*100; % cm
%             nden=P(i)/kb/T(i)*1E-6; % cm-3
%             clouddepth=clouddepth+(cloudopacity(i)+cloudmopacity(i))*log(P(i+1)/P(i))*sh; %cloud
%             gasdepth=gasdepth+(crossCH4*fCH4+crossH2O*watermix(i)+crossNH3*ammoniamix(i))*nden*log(P(i+1)/P(i))*sh;
%             rayleighdepth=rayleighdepth+crossRay*nden*log(P(i+1)/P(i))*sh;
%             i=i+1;
%         end
%         disp(P(i-1));
%         omega=(clouddepth+rayleighdepth)./(clouddepth+rayleighdepth+gasdepth);
%         albedoD=omega./(1+sqrt(1-omega)).^2;
%         
%         % Consider Partial Cloud Coverage
%         % methanedepth=P(end)/g*1.0E-4/MMM/AMU*crossCH4*fCH4;
%         % rayleighdepth=P(end)/g*1.0E-4/MMM/AMU*crossRay;
%         omega=(rayleighdepth)./(rayleighdepth+gasdepth);
%         albedoD0=omega./(1+sqrt(1-omega)).^2;
        
        % combine albedo
        albedoD = (albedoD*cloudfrac+albedoD0*(1-cloudfrac))/(2/3)*0.55; % 17.5% reduction due to Raman/Haze
        albedo = albedoD*solar/sum(solar)*0.8; % the phase integral = 0.67+-0.06 from 4 SS giant planets % assumed to be 1.25 Marley et al. (1999)
        
        disp(albedo)
        
        %albedo = 0.343;
        
        % Update LoopVar
        LoopVar = max(abs(Tnew-T)./Tnew);
        T = Tnew;
        % Iterate
        LoopID = LoopID+1;
    end
    
    % Standard Output
    %mkdir(['Result/',DIRName]);
    outdir=[DIRName,'/'];
    
    % Calculate the height
    P=fliplr(P);
    T=fliplr(T);
    watermix=fliplr(watermix);
    ammoniamix=fliplr(ammoniamix);
    sulfidemix=fliplr(sulfidemix);
    cloudden=fliplr(cloudden);
    cloudmden=fliplr(cloudmden);
    particlesize=fliplr(particlesize);
    particlemsize=fliplr(particlemsize);
    cloudcompo=fliplr(cloudcompo);
    Z(1)=0;
    for j=1:length(P)-1
        H = kb*(T(j)+T(j+1))/2/g/MMM/AMU/1000; %km
        Z(j+1) = Z(j) + H*log(P(j)/P(j+1));
    end
    
    % find tropopause
    tempgrad=abs(diff(log(T)));
    tempgradmax=max(tempgrad);
    i=length(Z);
    while tempgrad(i-1)<tempgradmax*0.3
        i = i-1;
    end
    idtropopause=i;
    ptropopause=P(idtropopause);
    while P(i)<P(idtropopause)*10
        i = i-1;
    end
    idbottom=i;
    % disp([P(idbottom),P(idtropopause)]);

    % Re-sample in equal height spacing
    nn = ceil(70/(Z(end)-Z(idbottom))*Z(end));
    zz = linspace(Z(1),Z(end),nn);
    T  = interp1(Z,T,zz);
    P  = exp(interp1(Z,log(P),zz));
    z0 = zz(1:nn-1);
    z1 = zz(2:nn);
    zl = z0*0.5+z1*0.5;
    tl = T(1:nn-1)*0.5+T(2:nn)*0.5;
    pl = sqrt(P(1:nn-1).*P(2:nn));
    nden=pl./kb./tl*1.0E-6; % molecule cm-3
    nH2 = fH2*nden;
    nHe = fHe*nden;
    nCH4= fCH4*nden;
    nH2O= exp(interp1(Z,log(watermix),zl)).*nden;
    nNH3= exp(interp1(Z,log(ammoniamix),zl)).*nden;
    nH2S= exp(interp1(Z,log(sulfidemix),zl)).*nden;
    cloudden = exp(interp1(Z,log(cloudden),zl));
    cloudmden = exp(interp1(Z,log(cloudmden),zl));
    particlesize = exp(interp1(Z,log(particlesize),zl));
    particlemsize = exp(interp1(Z,log(particlemsize),zl));
    cloudcompo = exp(interp1(Z,log(cloudcompo),zl));
    
    loglog(nH2O./nden,pl/1e+5,'k',nNH3./nden,pl/1e+5,'r',nH2S./nden,pl/1e+5,'c','LineWidth',2);
    hold on
    loglog(cloudden,pl/1e+5,'b',cloudmden,pl/1e+5,'m','LineWidth',2);
    loglog(tl/1000,pl/1e+5,'g');

    %loglog(exp(interp1(Z,log(watermix),zl)),pl,exp(interp1(Z,log(ammoniamix),zl)),pl);

    idbottom = nn-70;
    zbottom = zz(idbottom);
    idtropopause = round(interp1(P,1:1:nn,ptropopause));
    
    % Generate TP file
    dlmwrite([outdir,'TP.dat'],[zz(idbottom:end)-zbottom;log10(P(idbottom:end));T(idbottom:end)]','delimiter',' ','precision','%.6f');
    
    % Generate Kzz file
    % kzzt = (SB*Tint^4*kb./(P/kb./T*MMM*AMU)/MMM/AMU/14e+3).^(1/3).*(kb*T/MMM/AMU/g)*1E+4; % cm2/s
    kzzt = 1E+3*P./P; % cm2/s
    kzz_tropopause = kzzt(idtropopause);
    for i=1:length(zz)
        if i<=idtropopause
            kzz(i)=kzzt(i);
        else
            kzz(i)=kzzt(idtropopause)*(P(i)/T(i)*T(idtropopause)/P(idtropopause))^(-0.5);
        end
    end
    dlmwrite([outdir,'KZZ.dat'],[zz(idbottom:end)-zbottom;kzz(idbottom:end)]','delimiter',' ','precision','%.6f');

    % Generate ConcentrationSTD.dat file
    NSP=111;
    f = fopen([outdir,'ConcentrationSTD_R.dat'],'w');
    fprintf(f,'%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t','z','z0','z1','T','P');
    for i=1:NSP
        fprintf(f,'%d\t\t',i);
    end
    fprintf(f,'\n');
    fprintf(f,'%s\t\t%s\t\t%s\t\t%s\t\t%s\n','km','km','km','K','Pa');
    for j=1:length(zl)
        fprintf(f,'%.6f\t%.6f\t%.6f\t%.6f\t%.6e\t',zl(j),z0(j),z1(j),tl(j),pl(j));
        for i=1:NSP
            if i==7
                fprintf(f,'%.6e\t',nH2O(j));
            elseif i==21
                fprintf(f,'%.6e\t',nCH4(j));
            elseif i==9
                fprintf(f,'%.6e\t',nNH3(j));
            elseif i==45
                fprintf(f,'%.6e\t',nH2S(j));
            elseif i==53
                fprintf(f,'%.6e\t',nH2(j));
            else
                fprintf(f,'%.6e\t',0);
            end
        end
        fprintf(f,'\n');
    end
    fclose(f);
    
    f = fopen([outdir,'ConcentrationSTD_C.dat'],'w');
    fprintf(f,'%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t','z','z0','z1','T','P');
    for i=1:NSP
        fprintf(f,'%d\t\t',i);
    end
    fprintf(f,'\n');
    fprintf(f,'%s\t\t%s\t\t%s\t\t%s\t\t%s\n','km','km','km','K','Pa');
    for j=idbottom:length(zl)
        fprintf(f,'%.6f\t%.6f\t%.6f\t%.6f\t%.6e\t',zl(j)-zbottom,z0(j)-zbottom,z1(j)-zbottom,tl(j),pl(j));
        for i=1:NSP
            if i==7
                fprintf(f,'%.6e\t',nH2O(j));
            elseif i==21
                fprintf(f,'%.6e\t',nCH4(j));
            elseif i==9
                fprintf(f,'%.6e\t',nNH3(j));
            elseif i==45
                fprintf(f,'%.6e\t',nH2S(j));
            elseif i==53
                fprintf(f,'%.6e\t',nH2(j));
            else
                fprintf(f,'%.6e\t',0);
            end
        end
        fprintf(f,'\n');
    end
    fclose(f);
    
    % cloud output
    data = importdata('CrossP/cross_water_wavelength2.dat');
    H2OL_r = data(:,1); % zero-order radius, in micron
    H2OL_c = data(:,2:1388); % cross section per droplet, in cm2
    data = importdata('CrossP/albedo_water_wavelength2.dat');
    H2OL_a = data(:,2:1388);
    data = importdata('CrossP/geo_water_wavelength2.dat');
    H2OL_g = data(:,2:1388);
    
    data = importdata('CrossP/cross_water_wavelength2.dat');
    H2OI_r = data(:,1); % zero-order radius, in micron
    H2OI_c = data(:,2:1388); % cross section per droplet, in cm2
    data = importdata('CrossP/albedo_water_wavelength2.dat');
    H2OI_a = data(:,2:1388);
    data = importdata('CrossP/geo_water_wavelength2.dat');
    H2OI_g = data(:,2:1388);
    
    data = importdata('CrossP/cross_ammonia_wavelength2.dat');
    NH3I_r = data(:,1);  
    NH3I_c = data(:,2:1388);
    data = importdata('CrossP/albedo_ammonia_wavelength2.dat');
    NH3I_a = data(:,2:1388);
    data = importdata('CrossP/geo_ammonia_wavelength2.dat');
    NH3I_g = data(:,2:1388);
    
    % opacity
    sig=2;
    for j=1:length(zl)
        r2 = particlemsize(j);
        if cloudmden(j)<1e-20
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
        if cloudden(j)<1e-20
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
    
    f = fopen([outdir,'Cloudtop.dat'],'w');
    %Prefref=interp1(wave,Pref,800);
    fprintf(f,'Cloud top Pressure is at %.2e\n', Pref);
    fprintf(f,'Water cloud top Pressure is at %.2e\n', wPref);
    
    % add haze opacity
    hazetau=zeros(1,1387);
    for j=1:length(zl)
        if pl(j)>=Pref*exp(-2) && pl(j)<=Pref
            hazeden=hazemix*0.017*pl(j)/RGAS/tl(j); % kg/m^3, g/L
            [r0,r1,r2,VP]=particlesizef(g,tl(j),pl(j),MMM,17.0,KE,0.01*hazemix*pl(j));
            hazeopa=hazeden/VP*1.0E-3*10.^interp1(log10(NH3I_r),log10(NH3I_c),log10(max(0.01,min(r0,100)))); % cm-1
            hazealb=interp1(log10(NH3I_r),NH3I_a,log10(max(0.01,min(r0,100))));
            hazegeo=interp1(log10(NH3I_r),NH3I_g,log10(max(0.01,min(r0,100))));
            hazetau=hazetau+hazeopa*(z1(j)-z0(j))*1E+5;
            
            if hazemix>0
            alba(j,:)=(alba(j,:).*croa(j,:)+hazealb.*hazeopa)./(croa(j,:)+hazeopa);
            geoa(j,:)=(geoa(j,:).*croa(j,:)+hazegeo.*hazeopa)./(croa(j,:)+hazeopa);
            croa(j,:)=croa(j,:)+hazeopa;
            end
            
        end
    end
    
    lam1 = 0.1;
    while lam1(end)<100
        lam1(end+1)=lam1(end)*(1+1/200);
    end
   
    fprintf(f,'Photochemical haze has optical depth of %.2e\n', interp1(lam1,hazetau,0.8));
    fclose(f);
    
    dlmwrite([outdir,'Cross_H2O_R.dat'],crow,'delimiter','\t','precision','%.6e');
    dlmwrite([outdir,'albedo_H2O_R.dat'],albw,'delimiter','\t','precision','%.6e');
    dlmwrite([outdir,'geo_H2O_R.dat'],geow,'delimiter','\t','precision','%.6e');
    dlmwrite([outdir,'Cross_NH3_R.dat'],croa,'delimiter','\t','precision','%.6e');
    dlmwrite([outdir,'albedo_NH3_R.dat'],alba,'delimiter','\t','precision','%.6e');
    dlmwrite([outdir,'geo_NH3_R.dat'],geoa,'delimiter','\t','precision','%.6e');
    
    dlmwrite([outdir,'Cross_H2O_C.dat'],crow(idbottom:end,:),'delimiter','\t','precision','%.6e');
    dlmwrite([outdir,'albedo_H2O_C.dat'],albw(idbottom:end,:),'delimiter','\t','precision','%.6e');
    dlmwrite([outdir,'geo_H2O_C.dat'],geow(idbottom:end,:),'delimiter','\t','precision','%.6e');
    dlmwrite([outdir,'Cross_NH3_C.dat'],croa(idbottom:end,:),'delimiter','\t','precision','%.6e');
    dlmwrite([outdir,'albedo_NH3_C.dat'],alba(idbottom:end,:),'delimiter','\t','precision','%.6e');
    dlmwrite([outdir,'geo_NH3_C.dat'],geoa(idbottom:end,:),'delimiter','\t','precision','%.6e');
    
%     f = fopen([outdir,'Cross_H2O_C.dat'],'w');
%     for i=idbottom:length(zl)
%         for j=1:1387
%             fprintf(f,'%.6e\t',crow(i,j));
%         end
%         fprintf(f,'\n');
%     end
%     fclose(f);
%     f = fopen([outdir,'Cross_NH3_C.dat'],'w');
%     for i=idbottom:length(zl)
%         for j=1:1387
%             fprintf(f,'%.6e\t',croa(i,j));
%         end
%         fprintf(f,'\n');
%     end
%     fclose(f);
%     
%     f = fopen([outdir,'albedo_H2O_C.dat'],'w');
%     for i=idbottom:length(zl)
%         for j=1:1387
%             fprintf(f,'%.6f\t',albw(i,j));
%         end
%         fprintf(f,'\n');
%     end
%     fclose(f);
%     f = fopen([outdir,'albedo_NH3_C.dat'],'w');
%     for i=idbottom:length(zl)
%         for j=1:1387
%             fprintf(f,'%.6f\t',alba(i,j));
%         end
%         fprintf(f,'\n');
%     end
%     fclose(f);
%     
%     f = fopen([outdir,'geo_H2O_C.dat'],'w');
%     for i=idbottom:length(zl)
%         for j=1:1387
%             fprintf(f,'%.6f\t',geow(i,j));
%         end
%         fprintf(f,'\n');
%     end
%     fclose(f);
%     f = fopen([outdir,'geo_NH3_C.dat'],'w');
%     for i=idbottom:length(zl)
%         for j=1:1387
%             fprintf(f,'%.6f\t',geoa(i,j));
%         end
%         fprintf(f,'\n');
%     end
%     fclose(f);
    
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

end

