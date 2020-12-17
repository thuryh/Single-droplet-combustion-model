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
    MW_N = 28.01e-3; % kg/mol N2
    MW_O = 32e-3; %kg/mol O2
    MW_CO2 = 44e-3; %kg/mol CO2
    MW_H2O = 18e-3; %kg/mol H2O
    
    T_boil_F1 = 139 + 273.15; % K boiling point of m-Xylene
    T_boil_F2 = 228 + 273.15; % K boiling point of 2-EHA    
 
    T_inf = 300; % temperature at the infinity
    cp = (2*844+3*1930+3*4.76*1040)/9.76*5; % J/kg/K gas heat capacity, initial guess
    lambda_g = 0.2;  % W/m/K heat conductivity of air gases, initial gas
    % gas component analysis
    % C H O combustible
    Chemformu1 = [8,10,0];
    Chemformu2 = [8,16,2];
    % Solvent 1 Solvent + ()O2 -> ST11 CO2 + ST12 H2O;
    ST11 = Chemformu1(1);
    ST12 = Chemformu1(2)/2;
    % Solute 1 Solute + () O2 -> ST21 CO2 + ST22 H2O;
    ST21 = Chemformu2(1);
    ST22 = Chemformu2(2)/2;
    
    HvH2O = 40.65e3;%latent heat for H2O, 25oC, J/mol
    
    Qv1 = 42.65e3/MW_F1;% latent energy for m-Xylene evaporation J/kg
    Qv2 = 75.60e3/MW_F2;% latent energy for 2-EHA evaporation J/kg
    Qc1 = 4.3745e6/MW_F1 + 1/MW_F1*ST12*HvH2O; % standard combustion enthaphy for m-Xylene J/kg
    Qc2 = 4.5233e6/MW_F2 + 1/MW_F2*ST22*HvH2O; % standard combustion enthaphy for C8H16O2 J/kg; meanless for Ni(OH)2
    
    nu1 = (Chemformu1(1) + Chemformu1(2)/4 - Chemformu1(3)/2)*MW_O/MW_F1; % for species 1, m-Xylene
    nu2 = (Chemformu2(1) + Chemformu2(2)/4 - Chemformu2(3)/2)*MW_O/MW_F2; % for species 2, 2-EHA
   
    %% liquid phase
    lambda_l1 = 0.1032; % heat conductivity of m-Xylene, W/m/K
    lambda_l2 = 200.7e-3;% heat conductivity of 2EHA, W/m/K
    cp_l1 = 184.5/MW_F1; % liquid capacity of m-Xylene, J/kg/K
    cp_l2 = 309.75/MW_F2*2; % liquid capacity of 2-EHA, J/kg/K
    rho_l1 = 0.86*1e3; % liquid density of m-xylene, kg/m3
    rho_l2 = 903; % liquid density of 2-EHA, kg/m3
    MW1 = MW_F1; % m-xylene 
    MW2 = MW_F2; % 2-EHA

    Le = 10;
    Le_g = 1;