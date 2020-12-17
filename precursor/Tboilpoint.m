function [Xb] = Tboilpoint(Y,phi,rho_d,T_boil_F1,T_boil_F2,T_boil_Pr,T)
    global Qv1 Qv2 Qv3 MW_F1 MW_F2 MW_Pr Rg rho_s MW_S
    global X1 X2 X3 beta1 beta2 beta3
    global N
    Y1 = Y(1:N+1);
    Y2 = Y(N+2:2*N+2);
    Y3 = Y(2*N+3:3*N+3);
    
    beta1 = Qv1*MW_F1/Rg;
    beta2 = Qv2*MW_F2/Rg;
    beta3 = Qv3*MW_Pr/Rg;
    X1 = Y1./MW_F1./(Y1./MW_F1 + Y2./MW_F2 + Y3./MW_Pr);
    X2 = Y2./MW_F2./(Y1./MW_F1 + Y2./MW_F2 + Y3./MW_Pr);
    X3 = Y3./MW_Pr./(Y1./MW_F1 + Y2./MW_F2 + Y3./MW_Pr);
    Xb = X1.*exp(beta1*(1/T_boil_F1-1./T)) + X2.*exp(beta2*(1/T_boil_F2-1./T)) + X3.*exp(beta3*(1/T_boil_Pr-1./T));
end