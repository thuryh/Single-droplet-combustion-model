function Tb = Tboilpoint(Y,T_boil_F1,T_boil_F2)
    global Qv1 Qv2 MW_F1 MW_F2 Rg
    global X1 X2 beta1 beta2
    global N
    Y1 = Y(1:N+1);
    Y2 = Y(N+2:2*N+2);
    
    beta1 = Qv1*MW_F1/Rg;
    beta2 = Qv2*MW_F2/Rg;
    
    X1 = Y1./MW_F1./(Y1./MW_F1 + Y2./MW_F2);
    X2 = Y2./MW_F2./(Y1./MW_F1 + Y2./MW_F2);
    
    Tb0 = (X1.*beta1 + X2.*beta2)./(X1.*beta1/T_boil_F1 + X2.*beta2/T_boil_F2+X1+X2-1);
    options = optimoptions(@fmincon,'MaxFunctionEvaluations',1e5,'StepTolerance',1);
    [Tb,eval] = fsolve(@fun,Tb0,options);
    fprintf('error = %f\n',eval);
end
function y = fun(Tb)
global T_boil_F1 T_boil_F2 X1 X2 beta1 beta2 N
    y = sum(abs(X1.*exp(beta1*(1/T_boil_F1-1./Tb))+X2.*exp(beta2*(1/T_boil_F2-1./Tb))-1))/(N+1);
end
