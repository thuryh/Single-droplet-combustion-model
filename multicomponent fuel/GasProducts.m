function [MF_CO2,MF_H2O,XF_CO2,XF_H2O] = GasProducts(epsilon,ST11,ST12,ST21,ST22)
%Mass Fraction of CO2 H2O
global MW_F1 MW_F2
    %Mole of CO2 and H2O per molar reaction
    MW_CO2 = 44e-3;
    MW_H2O = 18e-3;
    CO2_Numper = epsilon*ST11/MW_F1 + (1 - epsilon)*ST21/MW_F2;
    H2O_Numper = epsilon*ST12/MW_F1 + (1 - epsilon)*ST22/MW_F2;
    
    MF_CO2 = CO2_Numper*MW_CO2/(CO2_Numper*MW_CO2 + H2O_Numper*MW_H2O);
    MF_H2O = H2O_Numper*MW_H2O/(CO2_Numper*MW_CO2 + H2O_Numper*MW_H2O);
    XF_CO2 = CO2_Numper/(CO2_Numper + H2O_Numper);
    XF_H2O = H2O_Numper/(CO2_Numper + H2O_Numper);
end

