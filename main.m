clear all; clc; close all;
global kb Rg Patm Na
global N xc xn dx r dr
global nu Y_O_inf MW_F MW_O MW_N MW_g T_boil_F Vp0 Qv Qc T_inf cp lambda_g rho_l rho_p K1 lambda_l dTdxs dYFdxs drs2dt alpha_l alpha Bm
global Le_l1 Le_l2 dY1dxs dY2dxs epsilon Ts Le_g MW_F1 MW_F2 MW_P

%% Grid
    N = 200;
    xc= linspace(0.5/N,1-0.5/N,N)';
    xn= linspace(0,1,N+1)';
    dx= 1/N;
    
    rmax = 50;
    rNum = 200;
    rNum = 1*rNum;
    r = linspace(1,rmax,rNum);%gas phase spatial resolution
    dr = (rmax-1)/(rNum-1);
    
    tspan = linspace(0,25e-3,1e3+1); 
%% volatle setting
    vola = 0;
    if vola == 0
       parametersetting_unvola;
    else
       parametersetting;
    end
    
    
%% preparing for initial conditions for gas phase equations variations
    Ts = 300;
    rs0 = 100e-6;%m
    rs = rs0;

%% numerical simulation setting
    
    T0 = ones( N+1,1)*Ts;
    Y0 = zeros(2*(N+1),1);
    Y0(N+2:2*N+2) = 23.792e-2;
    Y0(1:N+1) = 1 - Y0(N+2:2*N+2);
    phi0 = ones(N+1,1);
    ppsi0 = zeros(N+1,1);
    theta0 = zeros(N+1,1);
    
    Time_history=[];
    Tresult_history=[];
    Tgresult_history=[];
    Y1result_history=[];
    Y2result_history=[];
    phiresult_history = [];% liquid in Liquid phase
    ppsiresult_history = [];% precipitation in in Liquid phase
    thetaresult_history = [];% decompostion in Liquid phase
    YFresult_history=[];
    psiresult_history=[];
    nresult_history=[];
    vresult_history=[];
    dpresult_history = [];
    dmdt_history = [];
    rs_history = [];
    rf_history = [];
    drs2dt_history = [];
    B_history = [];
    Tinit = T0;
    Yinit = Y0;
    phiinit = phi0;
    ppsiinit = ppsi0;
    thetainit = theta0;
    T = ones(1,rNum)*300;% estimanted for gas phase
    Sp = zeros(1,rNum);
    j0 = 2;
    dTinfdt = 10^5;
    for i=1:1:length(tspan)-1
        for j = 1:j0
        %T_inf = T_inf+dTinfdt*(tspan(i+1)-tspan(i));
        %% liquid T and c for the gas-phase b.c. in the next step
        Ts  = Tinit(end,end); % droplet surface temperature
        Y1s = Yinit(1*N+1,end);% droplet     volatile species, mass fraction
        Y2s = Yinit(2*N+2,end);% droplet non-volatile species, mass fraction

        X1s = Y1s/MW1 / (Y1s/MW1 + Y2s/MW2); % droplet     volatile species, mole fraction
        X2s = Y2s/MW2 / (Y1s/MW1 + Y2s/MW2); % droplet non-volatile species, mole fraction
        
        XF1s = X1s*exp(Qv1*MW_F1/Rg*(1/T_boil_F1-1/Ts)); % volatile species,        mole fraction at the droplet surfac
        XF2s = X2s*exp(Qv2*MW_F2/Rg*(1/T_boil_F2-1/Ts));

        XNs = 1 - XF1s - XF2s;                                 % non-combustible species, mole fraction at the droplet surface
        YF1s = XF1s*MW_F1/(XF1s*MW_F1 + XF2s*MW_F2 + XNs*MW_N);
        YF2s = XF2s*MW_F2/(XF1s*MW_F1 + XF2s*MW_F2 + XNs*MW_N);
        YFs = YF1s + YF2s;
        if YFs > 1
            YF1s = YF1s/YFs;
            YF2s = YF2s/YFs;
            YFs = YF1s + YF2s;
        end

        epsilon = YF1s/YFs;
        MW_P = ((epsilon/MW_F1*ST11 + (1 - epsilon)/MW_F2*ST21)*MW_CO2 + (epsilon/MW_F1*ST12 + (1 - epsilon)/MW_F2*ST22)*MW_H2O)/(epsilon/MW_F1*(ST11 + ST12) + (1-epsilon)/MW_F2*(ST21 + ST22));
        %% mass and heat diffusion coefficient
        
        % averaged density in droplet
        rho_l = 1./(Yinit(1:N+1)/rho_l1 + Yinit(N+2:2*N+2)/rho_l2);
        rho_d = rho_l;
            % rho_d = 1./(phiinit./rho_l + ppsiinit./rho_pr + thetainit./rho_de);
            % heat diffusion cofficient
        cp_l = cp_l1*Yinit(1:N+1) + cp_l2*Yinit(N+2:2*N+2);
        
            % cp_d = cp_l.*phiinit + cp_pr.*ppsiinit + cp_de.*thetainit;
        cp_d = cp_l;
        lambda_l = 1./(Yinit(1:N+1)/lambda_l1 + Yinit(N+2:2*N+2)/lambda_l2);
        lambda_d = lambda_l;
            % lambda_d = 1./(phiinit./lambda_l + ppsiinit./lambda_pr + thetainit./lambda_de);
        alpha_d = lambda_d./cp_d./rho_d;
        
        D = alpha_d./Le;
        D = [D,D];
        
        %% gas phase temperature and species profile
        nu = epsilon*nu1 + (1 - epsilon)*nu2;
        Qc = epsilon*Qc1 + (1 - epsilon)*Qc2;
        Qv = epsilon*Qv1 + (1 - epsilon)*Qv2;
        Bm = (Y_O_inf/nu+YFs)/(1-YFs);
 
        CombCri = rs^2/rs0^2*100;
        CombCri0 = 0;
        if (CombCri < CombCri0)
            rf = 1;
        else
            rf = log(1+Bm)/log(1 +Y_O_inf/nu);
        end
        
        
        if (CombCri < CombCri0)
            YF = zeros(1,length(r));
            YO = zeros(1,length(r));
            YP = ones(1,length(r));
            YN = zeros(1,length(r));
        else
            YF(find(r< rf)) = -Y_O_inf/nu + (YFs+Y_O_inf/nu)*(exp(-1./r(find(r<rf))*log(1+Bm))-1)/(1/(1+Bm)-1);%1 - (Y+Y_O_inf/nu)*exp(-log(1+Bm)./r(find(r< rf)));
            YF(find(r>=rf)) = zeros(1,length(find(r>=rf)));
            YO(find(r< rf)) = zeros(1,length(find(r< rf)));
            YO(find(r>=rf)) = Y_O_inf - (YFs*nu+Y_O_inf)*(exp(-1./r(find(r>=rf))*log(1+Bm))-1)/(1/(1+Bm)-1);
            YP(find(r< rf)) = (nu+1)/nu*Y_O_inf - (nu+1)/nu*Bm/(1+Bm)*Y_O_inf*(exp(-1./r(find(r< rf))*log(1+Bm))-1)/(1/(1+Bm)-1);
            YP(find(r>=rf)) = ((nu+1)*YFs+(1+nu)/nu/(1+Bm)*Y_O_inf)*(exp(-1./r(find(r>=rf))*log(1+Bm))-1)/(1/(1+Bm)-1);
            YN              = Y_N_inf - Bm/(1+Bm)*Y_N_inf*(exp(-1./r*log(1+Bm))-1)/(1/(1+Bm)-1);
        end
        
        %% Physical parameter estimation
        rfp = ceil((rf - 1)/dr) + 1;% find flame position in grid
        [lambda_g,cp] = Gasmixpro(epsilon,ST11,ST12,ST21,ST22,YF,YO,YP,YN,T,rfp);

        T(find(r< rf)) = T_inf + Y_O_inf*Qc/nu/cp + (Ts - T_inf - Y_O_inf*Qc/nu/cp)*(exp(-1./r(find(r< rf))*log(1+Bm))-1)/(1/(1+Bm)-1);
        T(find(r>=rf)) = T_inf + (Ts - T_inf + YFs*Qc/cp)*(exp(-1./r(find(r>=rf))*log(1+Bm))-1)/(1/(1+Bm)-1);

        %% Parameter
        drs2dt = -2*lambda_g/cp/rho_d(end)*log(1+Bm)/Le_g;

        %% liquid prop
        Tresult = heattransfer(linspace(tspan(i)*(j0+1-j)/j0 + tspan(i+1)*(j-1)/j0,tspan(i+1),3),Tinit,rho_d,rs,alpha_d,cp_d,drs2dt);
        Yresult = masstransfer(linspace(tspan(i)*(j0+1-j)/j0 + tspan(i+1)*(j-1)/j0,tspan(i+1),3),Yinit,rho_d,Tinit,rs,D,drs2dt);

        Y1result  = Yresult(1:N+1);
        Y2result  = 1 - Y1result;
        
        if j ~= j0
            jj = (j0-j)/(j0+1-j);
            Tinit = jj.*Tresult + (1 - jj).*Tinit;
            Yinit = jj.*[Y1result;Y2result] + (1 - jj).*Yinit;
        else
            Tinit = Tresult;
            Yinit = [Y1result;Y2result];
        end
        
        end
        YFresult_history=[YFresult_history;[YF1s,YF2s,YF1s+YF2s]];
        
        rf_history = [rf_history;rf];
        
        Tgresult_history = [Tgresult_history;T];
        
        Tresult_history   = [Tresult_history;Tresult'];
        Y1result_history  = [Y1result_history;Y1result'];
        Y2result_history  = [Y2result_history;Y2result'];

        
        
        
        drs2dt_history = [drs2dt_history;drs2dt];
        B_history = [B_history;Bm];
        
        %% solve droplet shrink
        dmdt = 4*pi*rs*lambda_g/cp*log(1+Bm)/Le_g; % changing mass flow rate, kg/s

        
        rs2 = rs^2+drs2dt*(tspan(i+1)-tspan(i));
        if ((1 - Y1s)*(Y1s - 0) < 0 )
            break
        else
            if((1 - Y2s)*(Y2s - 0) < 0)
                break
            else
                if ((1 - YF1s)*(YF1s - 0) < 0)
                    break
                else
                    if ((1 - YF2s)*(YF2s - 0) < 0)
                        break
                    end
                end
            end
        end
        

        if (rs2 > 0)
            rs = sqrt(rs2);
            rs_history = [rs_history;rs];
        else
            break
        end
        if (CombCri < CombCri0)
            break
        end
      
        %% solve the PBM model
        MW_g = (epsilon*YF/MW_F1 + (1 - epsilon)*YF/MW_F2 + YO/MW_O + YP/MW_P + YN/MW_N).^(-1);% kg mol-1
        rho_g = Patm.*MW_g/Rg./T; % kg m-3
        dmdt = 4*pi*rs*lambda_g/cp*log(1+Bm)/Le_g; % changing mass flow rate, kg/s
        dmdt_history = [dmdt_history;dmdt];
        thermal_velo;
        v = dmdt./rho_g./(4*pi*r.^2*rs^2) + vth;% gas velocity, m/s
        vresult_history = [vresult_history;v];
        rfposition = ceil((rf - 1)/dr) + 1;% find flame position in grid
        Sp(rfposition) = 0*(dmdt*YF2s - (1 - epsilon)*4*pi*rs^2*lambda_g/cp/Le_g*(YF2s + Y_O_inf/nu)*(-log(1 + Bm)/Bm/rs))/(4*pi*rf^2*rs^2)*M_mon/(M_pre*rho_p)*nu_pre/dr/rs;% source term m s-1
        PBM = [psi,n];
        [Time,PBMresult] = ode45(@(t,PBM)PBMsolve(t,PBM,Sp,rs,v,rho_g,T,drs2dt),linspace(tspan(i),tspan(i+1),3),PBM);
        
        psi = PBMresult(end,1:rNum);
        n   = PBMresult(end,rNum+1:end);
        [vp,dp] = calvpdp(psi,n);

        psiresult_history = [psiresult_history;psi];
        nresult_history = [nresult_history;n];
        dpresult_history = [dpresult_history;dp];
        
        Sp(rfposition) = 0;
        
            
        Time_history=[Time_history;tspan(i+1)];
        %% temporal output 
        if(mod(i,10)==0)
        fprintf('step = %d, r^2 = %g%% \n',i,rs^2/rs0^2*100);
        figure(1);set(gcf, 'position', [500,0,1000,1000]);
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
        plot(Time_history,rs_history.^2/rs0^2)
        %yyaxis left; plot(r,psiresult_history(end,:),'b-');xlabel('Dimensionless radial position r/r_s');ylabel('Particle volume fraction');
        %yyaxis right;plot(r,dpresult_history(end,:)*1e9,'r-');xlabel('Dimensionless radial position r/r_s');ylabel('Particle size (nm)');
        set(gca,'FontSize',14);
        end


        

    end
