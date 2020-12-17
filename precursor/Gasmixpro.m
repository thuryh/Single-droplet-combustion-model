function [lambda_g,lambda_bd,cp,cp_bd,MW_P] = Gasmixpro(epsilon1,epsilon2,epsilon3,ST11,ST12,ST21,ST22,ST31,ST32,YF,YO,YP,YN,T,rfp)
    global MW_F1 MW_F2 Rg MW_P MW_Pr
    MW_N = 28.01e-3; % kg/mol, N2
    MW_O = 32e-3; %kg/mol, O2
    MW_CO2 = 44e-3; %kg/mol, CO2
    MW_H2O = 18e-3; %kg/mol, H2O
    [MF_CO2,MF_H2O,XF_CO2,XF_H2O] = GasProducts(epsilon1,epsilon2,epsilon3,ST11,ST12,ST21,ST22,ST31,ST32);

    YF1 = YF*epsilon1;
    YF2 = YF*epsilon2;
    YFPr = YF*epsilon3;
    YCO2= YP*MF_CO2;
    YH2O= YP*MF_H2O;

    XF1  = YF1/MW_F1  ./(YF1/MW_F1+YF2/MW_F2+YFPr/MW_Pr+YO/MW_O+YN/MW_N+YCO2/MW_CO2+YH2O/MW_H2O);
    XF2  = YF2/MW_F2  ./(YF1/MW_F1+YF2/MW_F2+YFPr/MW_Pr+YO/MW_O+YN/MW_N+YCO2/MW_CO2+YH2O/MW_H2O);
    XPr  = YFPr/MW_Pr ./(YF1/MW_F1+YF2/MW_F2+YFPr/MW_Pr+YO/MW_O+YN/MW_N+YCO2/MW_CO2+YH2O/MW_H2O);
    XO   = YO/MW_O    ./(YF1/MW_F1+YF2/MW_F2+YFPr/MW_Pr+YO/MW_O+YN/MW_N+YCO2/MW_CO2+YH2O/MW_H2O);
    XN   = YN/MW_N    ./(YF1/MW_F1+YF2/MW_F2+YFPr/MW_Pr+YO/MW_O+YN/MW_N+YCO2/MW_CO2+YH2O/MW_H2O);
    XCO2 = YCO2/MW_CO2./(YF1/MW_F1+YF2/MW_F2+YFPr/MW_Pr+YO/MW_O+YN/MW_N+YCO2/MW_CO2+YH2O/MW_H2O);
    XH2O = YH2O/MW_H2O./(YF1/MW_F1+YF2/MW_F2+YFPr/MW_Pr+YO/MW_O+YN/MW_N+YCO2/MW_CO2+YH2O/MW_H2O);
    XP   = XCO2 + XH2O;
    MW_P = (XCO2(1)*MW_CO2+XH2O(1)*MW_H2O)/(XCO2(1)+XH2O(1));
    % characteristic temperature and componient
    ratio1 = 1/6; ratio2 = 2/3; ratio3 = 1/6;
    Tch   = ratio1*T(1)   + ratio2*T(rfp)   + ratio3*T(end);
    YF1ch = ratio1*YF1(1) + ratio2*YF1(rfp) + ratio3*YF1(end);
    YF2ch = ratio1*YF2(1) + ratio2*YF2(rfp) + ratio3*YF2(end);
    YPrch = ratio1*YFPr(1)+ ratio2*YFPr(rfp)+ ratio3*YFPr(end);
    YOch  = ratio1*YO(1)  + ratio2*YO(rfp)  + ratio3*YO(end);
    YPch  = ratio1*YP(1)  + ratio2*YP(rfp)  + ratio3*YP(end);
    YNch  = ratio1*YN(1)  + ratio2*YN(rfp)  + ratio3*YN(end);
    YCO2ch = YPch*MF_CO2;
    YH2Och = YPch*MF_H2O;
    
    XF1ch = ratio1*XF1(1) + ratio2*XF1(rfp)+ ratio3*XF1(end);
    XF2ch = ratio1*XF2(1) + ratio2*XF2(rfp)+ ratio3*XF2(end);
    XPrch = ratio1*XPr(1) + ratio2*XPr(rfp)+ ratio3*XPr(end);
    XOch  = ratio1*XO(1)  + ratio2*XO(rfp) + ratio3*XO(end);
    XPch  = ratio1*XP(1)  + ratio2*XP(rfp) + ratio3*XP(end);
    XNch  = ratio1*XN(1)  + ratio2*XN(rfp) + ratio3*XN(end);
    XCO2ch = XPch*XF_CO2;
    XH2Och = XPch*XF_H2O;
    
%% Expression of thermal conductivity for specific heat and kinetic viscouty of gas species
    %% thermal conductivity
    Ych = [YF1ch,YF2ch,YPrch,YCO2ch,YOch,YNch,YH2Och];
    Xch = [XF1ch,XF2ch,XPrch,XCO2ch,XOch,XNch,XH2Och];
    lambda_g = lamb(Ych,Xch,Tch);
    Ybd = [YF1(1),YF2(1),YFPr(1),YCO2(1),YO(1),YN(1),YH2O(1)];
    Xbd = [XF1(1),XF2(1),XPr(1), XCO2(1),XO(1),XN(1),XH2O(1)];
    lambda_bd = lamb(Ybd,Xbd,T(1));

    %% specific heat J/kg/K
    cp_F1  = specheat('C8H10',Tch)/MW_F1;
    cp_F2  = specheat('C8H16O2',Tch)/MW_F2;
    cp_Pr  = specheat('C8H16O2',Tch)/MW_Pr; % the heat capacity of the precursor vapor is estimated by 2-EHA
    cp_O2  = specheat('O2' ,Tch)/MW_O;
    cp_CO2 = specheat('CO2',Tch)/MW_CO2;
    cp_H2O = specheat('H2O',Tch)/MW_H2O;
    cp_N2  = specheat('N2' ,Tch)/MW_N;
    cp = YF1ch*cp_F1 + YF2ch*cp_F2 + YPrch*cp_Pr+ YOch*cp_O2 + YNch*cp_N2 + YCO2ch*cp_CO2 + YH2Och*cp_H2O;
    
    cp_F1  = specheat('C8H10',T(1))/MW_F1;
    cp_F2  = specheat('C8H16O2',T(1))/MW_F2;
    cp_Pr  = specheat('C8H16O2',T(1))/MW_Pr; % the heat capacity of the precursor vapor is estimated by 2-EHA
    cp_O2  = specheat('O2' ,T(1))/MW_O;
    cp_CO2 = specheat('CO2',T(1))/MW_CO2;
    cp_H2O = specheat('H2O',T(1))/MW_H2O;
    cp_N2  = specheat('N2' ,T(1))/MW_N;
    cp_bd = YF1(1)*cp_F1 + YF2(1)*cp_F2 + YFPr(1)*cp_Pr + YO(1)*cp_O2 + YN(1)*cp_N2 + YCO2(1)*cp_CO2 + YH2O(1)*cp_H2O;

end