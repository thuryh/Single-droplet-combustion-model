clear all; clc; close all;
global kb Rg Patm Na
global N xc xn dx r dr
global nu Y_O_inf MW_F MW_O MW_N MW_g T_boil_F Vp0 Qv Qc T_inf cp lambda_g lambda_bd rho_l rho_p K1 lambda_l dTdxs dYFdxs drs2dt alpha_l alpha Bm
global Le_l1 Le_l2 dY1dxs dY2dxs epsilon Ts Le_g MW_F1 MW_F2 MW_P cp_bd Qv1 Qv2 T_boil_F1 T_boil_F2

%% Grid
N = 500; % droplet grid
xc= linspace(0.5/N,1-0.5/N,N)'; % central point in the droplet grid (non-dimensional)
xn= linspace(0,1,N+1)'; % node point in the droplet grid (non-dimensional)
dx= 1/N; % grid size in the droplet (non-dimensional)

rmax = 50; % the farthest position of the gas phase(non-dimensional)
rNum = 200;% the number of grid of the gas-phase flame
r = linspace(1,rmax,rNum); % grid of the gas phase (non-dimensional)
dr = (rmax-1)/(rNum-1); % grid size of the gas phase(non-dimensional)
tspan = linspace(0,25e-3,1e4+1); %grid of time (unit: s)

%% parameter setting
parametersetting_unvola;

%% preparing for initial conditions
Ts = 300; %surface temperature (K)
rs0 = 50e-6;%droplet size (m)
rs = rs0;%setting of temporal droplet size
rho_d = rho_l1;%initial droplet density is equal to the major fuel component
MW_P = MW_N;% molar mass for product species is equal to the inert species
T = ones(1,rNum)*300;% initial estimation temperature for gas phase (K)
XNs = 1;% mole fraction of inert species is initially set as 1
XPs = 0;% mole fraction of product species is initially set as 0

%% numerical simulation setting
T0 = ones( N+1,1)*Ts; %initial temperature
Y0 = zeros(2*(N+1),1);%initial mass fractions of components in droplets
Y0(N+2:2*N+2) = ones(N+1,1)*23.792e-2; % less-volatile component, 2-EHA
Y0(1:N+1) = 1 - Y0(N+2:2*N+2); % major fuel, m-xylene
Tinit = T0;
Yinit = Y0;

%% parameters for storage
Time_history=[];
Tresult_history=[];
Tgresult_history=[];
Y1result_history=[];
Y2result_history=[];
YFresult_history=[];
rs_history = [];
rf_history = [];
drs2dt_history = [];
B_history = [];

j0 = 5;
for i=1:1:length(tspan)-1
    for j = 1:j0
        %% liquid T and c for the gas-phase boundary condition
        Ts  = Tinit(end,end); % droplet surface temperature
        Y1s = Yinit(1*N+1,end);% droplet     volatile species, mass fraction
        Y2s = Yinit(2*N+2,end);% droplet non-volatile species, mass fraction

        X1s = Y1s/MW1 / (Y1s/MW1 + Y2s/MW2); % droplet     volatile species, mole fraction
        X2s = Y2s/MW2 / (Y1s/MW1 + Y2s/MW2); % droplet non-volatile species, mole fraction
        
        XF1s = X1s*exp(Qv1*MW_F1/Rg*(1/T_boil_F1-1/Ts)); % mole fraction of species 1 at the droplet surface in the gas phase
        XF2s = X2s*exp(Qv2*MW_F2/Rg*(1/T_boil_F2-1/Ts)); % mole fraction of species 2 at the droplet surface in the gas phase
        
        YF1s = XF1s*MW_F1/(XF1s*MW_F1 + XF2s*MW_F2 + XNs*MW_N + XPs*MW_P); % mass fraction of species 1 at the droplet surface in the gas phase
        YF2s = XF2s*MW_F2/(XF1s*MW_F1 + XF2s*MW_F2 + XNs*MW_N + XPs*MW_P); % mass fraction of species 2 at the droplet surface in the gas phase
        YFs = YF1s + YF2s; % mass fraction of the total fuel species at the droplet surface in the gas phase 
        epsilon = YF1s/YFs;% fractional mass vaporization rate of species 1
		
        MW_P = ((epsilon/MW_F1*ST11 + (1 - epsilon)/MW_F2*ST21)*MW_CO2 + (epsilon/MW_F1*ST12 + (1 - epsilon)/MW_F2*ST22)*MW_H2O)/(epsilon/MW_F1*(ST11 + ST12) + (1-epsilon)/MW_F2*(ST21 + ST22)); % molar mass of product species based on the evaporated fuel species (kg/mol)
        MW_F = 1/(epsilon/MW_F1+(1-epsilon)/MW_F2); % molar fraction of fuel species (kg/mol)
        
        %% gas phase temperature and species profile
        nu = epsilon*nu1 + (1 - epsilon)*nu2; % averaged stoichiometric oxygen-to-fuel mass ratio 
        Qc = epsilon*Qc1 + (1 - epsilon)*Qc2; % combustion heat (J/kg)
        Qv = epsilon*Qv1 + (1 - epsilon)*Qv2; % vaporization heat (J/kg)
        Bm = (Y_O_inf/nu+YFs)/(1-YFs); % Spalding transfer number 
 
        CombCri = rs^2/rs0^2*100; % The critical combustion limit, here we assume the flame exists throughout the droplet lifetime
        CombCri0 = 0;
        if (CombCri < CombCri0)
            rf = 1;
        else
            rf = log(1+Bm)/log(1 +Y_O_inf/nu); % flame front position, non-dimensional
            YF(find(r< rf)) = -Y_O_inf/nu + (YFs+Y_O_inf/nu)*(exp(-1./r(find(r<rf))*log(1+Bm))-1)/(1/(1+Bm)-1); %distribution of fuel species inside the flame front
            YF(find(r>=rf)) = zeros(1,length(find(r>=rf))); %distribution of fuel species outside the flame front
			YO(find(r< rf)) = zeros(1,length(find(r< rf))); %distribution of oxidizer inside the flame front
            YO(find(r>=rf)) = Y_O_inf - (YFs*nu+Y_O_inf)*(exp(-1./r(find(r>=rf))*log(1+Bm))-1)/(1/(1+Bm)-1); %distribution of oxidizer outside the flame front
			YP(find(r< rf)) = (nu+1)/nu*Y_O_inf - (nu+1)/nu*Bm/(1+Bm)*Y_O_inf*(exp(-1./r(find(r< rf))*log(1+Bm))-1)/(1/(1+Bm)-1);  %distribution of product inside the flame front
            YP(find(r>=rf)) = ((nu+1)*YFs+(1+nu)/nu/(1+Bm)*Y_O_inf)*(exp(-1./r(find(r>=rf))*log(1+Bm))-1)/(1/(1+Bm)-1); %distribution of product outside the flame front
			YN              = Y_N_inf - Bm/(1+Bm)*Y_N_inf*(exp(-1./r*log(1+Bm))-1)/(1/(1+Bm)-1); %distribution of the inert species
        end
        YFs = YF(1); % fuel mass fraction at the droplet surface
        YOs = YO(1); % oxidizer mass fraction at the droplet surface
        YPs = YP(1); % product species mass fraction at the droplet surface
        YNs = YN(1); % inert species mass fraction at the droplet surface
        XPs = YPs/MW_P/(YFs/MW_F + YOs/MW_O + YPs/MW_P + YNs/MW_N); % mole fraction of product species at the droplet surface
        XNs = YNs/MW_N/(YFs/MW_F + YOs/MW_O + YPs/MW_P + YNs/MW_N); % mole fraction of inert species at the droplet surface
        
        %% Physical parameter estimation
        rfp = ceil((rf - 1)/dr) + 1;% find flame position in grid
        [lambda_g,lambda_bd,cp,cp_bd,MW_P] = Gasmixpro(epsilon,ST11,ST12,ST21,ST22,YF,YO,YP,YN,T,rfp);% estimate the gas-phase thermal conductivity (J/m/K), heat capacity(J/kg/K), molar mass of product ('g' for average value, 'bd' for surface value), 
        T(find(r< rf)) = T_inf + Y_O_inf*Qc/nu/cp + (Ts - T_inf - Y_O_inf*Qc/nu/cp)*(exp(-1./r(find(r< rf))*log(1+Bm))-1)/(1/(1+Bm)-1); %distribution of temperature inside the flame front
        T(find(r>=rf)) = T_inf + (Ts - T_inf + YFs*Qc/cp)*(exp(-1./r(find(r>=rf))*log(1+Bm))-1)/(1/(1+Bm)-1); %distribution of temperature outside the flame front

        %% Parameter
        drs2dt = -2*lambda_g/cp/rho_d(end)*log(1+Bm); % dr^2/dt (m2/s)

        %% mass and heat diffusion coefficient in liquids
        converge = 1;
        Tresult_temp = Tinit;
        Yresult_temp = Yinit;
        while(converge == 1) % this iteration for the variation of thermoproperties of liquid species
            [rho_d, cp_d, lambda_d] = Liquidmixpro(Yresult_temp,Tresult_temp,cp_l1,cp_l2,lambda_l1,lambda_l2); %estimate the density (kg/m3), heat capacity (J/kg/K), thermal conductivity(J/m/K)
            alpha_d = lambda_d./cp_d./rho_d; % thermal diffusivity (m2/s)
            D = alpha_d./Le; % mass diffusivity (m2/s)
            D = [D,D];
            Tresult = heattransfer(linspace(tspan(i)*(j0+1-j)/j0 + tspan(i+1)*(j-1)/j0,tspan(i+1),3),Tinit,rs,alpha_d,lambda_d); %solve the temperature distribution in the droplet
            Yresult = masstransfer(linspace(tspan(i)*(j0+1-j)/j0 + tspan(i+1)*(j-1)/j0,tspan(i+1),3),Yinit,rho_d,rs,D);%solve the species mass fraction in the droplet
            Y1result  = Yresult(1:N+1); % mass fraction for species 1
            Y2result  = 1 - Y1result; % mass fraction for species 2
            errorT = sum(abs(Tresult - Tresult_temp))/sum(abs(Tresult_temp));
            errorY = sum(abs(Yresult - Yresult_temp))/sum(abs(Yresult_temp));
            fprintf('rs=%f,%e,%e\n',rs,errorT,errorY);
            if(errorT < 1e-7 && errorY < 1e-7)
                converge = 0;
                Y1result  = Yresult(1:N+1);
                Y2result  = Yresult(N+2:2*N+2);
            else
                Tresult_temp = Tresult;
                Yresult_temp = Yresult;
            end
        end
        Tinit = Tresult;
        Yinit = Yresult;
        if j ~= j0 % the relaxation coefficient is added here to avoid the fluctuations caused by the iteration between gas and liquid phase.
            jj = (j0-j)/(j0+1-j);
            Tinit = jj.*Tresult + (1 - jj).*Tinit;
            Yinit = jj.*[Y1result;Y2result] + (1 - jj).*Yinit;
        else
            Tinit = Tresult;
            Yinit = [Y1result;Y2result];         
        end
    end
	%% save the result to history
    YFresult_history=[YFresult_history;[YF1s,YF2s,YF1s+YF2s]];
    rf_history = [rf_history;rf];
    Tgresult_history = [Tgresult_history;T];
    Tresult_history   = [Tresult_history;Tresult'];
    Y1result_history  = [Y1result_history;Y1result'];
    Y2result_history  = [Y2result_history;Y2result'];
    drs2dt_history = [drs2dt_history;drs2dt];
    B_history = [B_history;Ts];	
	Time_history=[Time_history;tspan(i+1)];

    %% solve droplet shrink
    dmdt = 4*pi*rs*lambda_g/cp*log(1+Bm)/Le_g; % changing mass flow rate, kg/s
    rs2 = (rs^2+drs2dt*(tspan(i+1)-tspan(i)));%/alpha_droplet^(2/3);
    if (rs2 > 0)
        rs = sqrt(rs2);
        rs_history = [rs_history;rs];
    else
        break
    end
    if(rs2/rs0^2 < 0.01) % when the droplet is too small, stop simulation to avoid the error when the droplet radius is 0.
        break
    end  
    
    %% temporal output 
    if(mod(i,1)==0)
        fprintf('step = %d, r^2 = %g, Bm = %e \n',i,rs^2/rs0^2*100,Bm);
        figure(1);set(gcf, 'position', [500,-300,1000,1000]);
        subplot(2,2,1);
        yyaxis left; plot(r,T,'b-');           xlabel('Dimensionless radial position r/r_s');ylabel('Temperature T(K)');
        yyaxis right;plot(r,YO,'r-',r,YF,'r--');xlabel('Dimensionless radial position r/r_s');ylabel('Species mass fraction Y_i');
        set(gca, 'XScale', 'log','FontSize',14);title(['t =',num2str(tspan(i+1)),'s']);
        
        subplot(2,2,2);
        plot(r,YN,'b-',r,YP,'r-');xlabel('Dimensionless radial position r/r_s');ylabel('Species mass fraction');
        set(gca, 'XScale', 'log','FontSize',14);
        
        subplot(2,2,3);
        yyaxis left; plot(xn,Tresult_history(end,:),'b-');xlabel('Dimensionless radial position r/r_s');ylabel('Temperature T(K)');
        yyaxis right;plot(xn,Y1result_history(end,:),'m-',xn,Y2result_history(end,:),'g-');xlabel('Dimensionless radial position r/r_s');ylabel('Species mass fraction');
        set(gca, 'FontSize',14);

        subplot(2,2,4);
        yyaxis left;plot(Time_history/rs0^2/1e6/4,rs_history.^2/rs0^2)
        yyaxis right;plot(Time_history/rs0^2/1e6/4,B_history)
        set(gca,'FontSize',14);
    end
end