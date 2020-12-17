function [rho_d, cp_d, lambda_d] = Liquidmixpro(Yinit,Tinit,cp_l1,cp_l2,lambda_l1,lambda_l2)
    global N 
    % liquid density constant
    rho_l1 = 995.58-0.1449*Tinit-0.0011*Tinit.^2;%860;
    rho_l2 = 903;
    rho_d = 1./(Yinit(1:N+1)./rho_l1 + Yinit(N+2:2*N+2)./rho_l2);
    
    % liquid heat capacity
    cp_l1 = heatcapacityliquid('C8H10',Tinit);
    cp_l2 = heatcapacityliquid('C8H16O2',Tinit);
    cp_d = cp_l1.*Yinit(1:N+1) + cp_l2.*Yinit(N+2:2*N+2);
    
    % liquid thermal conductivity
    lambda_l1 = thermalconductivityliquid('C8H10', Tinit);
    lambda_l2 = thermalconductivityliquid('C8H16O2', Tinit);
    lambda_d = 1./(Yinit(1:N+1)./lambda_l1 + Yinit(N+2:2*N+2)./lambda_l2);
end
