global kb Rg Patm
global N xc xn dx
global nu Y_O_inf MW_F MW_O MW_N T_boil_F Qv Qc T_inf cp lambda_g rho_l K1 lambda_l dTdxs dYFdxs drs2dt
global dp0 vp0 rho_l2
global Le_l1 Le_l2 dY1dxs dY2dxs
    %% constant
    kb = 1.3806505e-23;
    Rg = 8.314;
    Patm = 101325;
    Na = 6.02e23;
    
    %% gas phase 
    Y_O_inf = 1; % infinity position for air
    Y_N_inf = 1-Y_O_inf;
    MW_F1 = 106.16e-3; % kg/mol m-Xylene 
    MW_F2 = 144.214e-3; %kg/mol 2-EHA
    MW_N = 28.01e-3; % kg/mol
    MW_O = 32e-3; %kg/mol
    MW_CO2 = 44e-3;
    MW_H2O = 18e-3;
    
    T_boil_F1 = 139 + 273.15; % K m-Xylene
    T_boil_F2 = 228 + 273.15; % K 2-EHA    
 
    T_inf = 300;
    cp = (2*844+3*1930+3*4.76*1040)/9.76*5; % J/kg/K
    lambda_g = 0.2;  % heat conductivity of air gases W/m/K
    % gas component analysis
    % C H O combustible
    Chemformu1 = [8,10,0];
    Chemformu2 = [8,16,2];
    % Solvent 1 Sovlent + ()O2 -> ST11 CO2 + ST12 H2O;
    ST11 = Chemformu1(1);
    ST12 = Chemformu1(2)/2;
    % Solute 1 Solute + () O2 -> ST21 CO2 + ST22 H2O;
    ST21 = Chemformu2(1);
    ST22 = Chemformu2(2)/2;
    
    HvH2O = 40.65e3;%latent heat for H2O, 25oC, J/mol
    
    Qv1 = 42.65e3/MW_F1;% latent energy for m-Xylene evaporation J/kg
    Qv2 = 75.60e3/MW_F2;% 5 times of latent energy for EHA evaporation; meanless for Ni(OH)2
    Qc1 = 4.3745e6/MW_F1 + 1/MW_F1*ST12*HvH2O; % standard combustion enthaphy for m-Xylene J/kg
    Qc2 = 4.5233e6/MW_F2 + 1/MW_F2*ST22*HvH2O; % standard combustion enthaphy for C8H16O2 J/kg; meanless for Ni(OH)2
    
    % Sovlent
    nu1 = (Chemformu1(1) + Chemformu1(2)/4 - Chemformu1(3)/2)*MW_O/MW_F1; % for C4H10O
    % Solute
    nu2 = (Chemformu2(1) + Chemformu2(2)/4 - Chemformu2(3)/2)*MW_O/MW_F2; % for C8H16O2; meanless for Ni(OH)2
   
    %% liquid phase
    lambda_l1 = 0.1032; % heat conductivity of m-Xylene   W/m/K
    lambda_l2 = 200.7e-3;% 2EHA, replaced by 2-Ethylhexanol 
    lambda_pr  = 0.3*lambda_l1;
    lambda_de  = 0.3*lambda_l1;%
    cp_l1 = 184.5/MW_F1; % liquid capacity J/kg/K
    cp_l2 = 309.75/MW_F2*2;
    cp_pr = cp_l2;
    cp_de = cp_l2;
    rho_l1 = 0.86*1e3; % liquid density kg/m3
    rho_l2 = 903;
    rho_pr = rho_l2;
    rho_de = rho_l2;
    MW1 = MW_F1; % fuel 
    MW2 = MW_F2; % vola
    MW3 = MW2;%precipate;
    MW4 = MW2;%decomp;

    Le = 10;
    Le_g = 1;
    
    %% PBModel
    dp0 = 0.3e-9;
    vp0 = 1/6*pi*dp0^3;
    dp = dp0*ones(1,rNum);%m
    Vp = pi*dp.^3/6;%m3
    delta = zeros(1,rNum);
    rho_p = 1000;
    alpha = 0.9;
    psi = zeros(1,rNum);
    n = zeros(1,rNum);
%     M_pre = MW_F2;% kg mol-1; molecular weight for precusor;
    M_pre =  MW_F2; %Precusor
    M_mon = 101.96e-3;% kg mol-1; molecular weight for monomers by precusor;
    nu_pre = 0.5;% stiochmetric number for precusor reaction PRE - > M? % Al2O3
%     Y_solubility = 1;% Ni(OH)2 in C2H5OH
    T_decomposition = 135 + 273.15;% decomposition temp for Ni(OH)2