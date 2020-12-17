function [MF_CO2,MF_H2O,XF_CO2,XF_H2O] = GasProducts(epsilon1,epsilon2,epsilon3,ST11,ST12,ST21,ST22,ST31,ST32)
%Mass Fraction of CO2 H2O
global MW_F1 MW_F2 MW_Pr
    %Mole of CO2 and H2O per molar reaction
    MW_CO2 = 44e-3;
    MW_H2O = 18e-3;
    CO2_Numper = epsilon1*ST11/MW_F1 + epsilon2*ST21/MW_F2 + epsilon3*ST31/MW_Pr;
    H2O_Numper = epsilon1*ST12/MW_F1 + epsilon2*ST22/MW_F2 + epsilon3*ST32/MW_Pr;
    
    MF_CO2 = CO2_Numper*MW_CO2/(CO2_Numper*MW_CO2 + H2O_Numper*MW_H2O);
    MF_H2O = H2O_Numper*MW_H2O/(CO2_Numper*MW_CO2 + H2O_Numper*MW_H2O);
    XF_CO2 = CO2_Numper/(CO2_Numper + H2O_Numper);
    XF_H2O = H2O_Numper/(CO2_Numper + H2O_Numper);
end

