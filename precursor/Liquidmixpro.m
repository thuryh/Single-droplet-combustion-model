function [rho_l, cp_l, lambda_l, viscosity_l] = Liquidmixpro(Yinit,Tinit)
    global N 
    alpha = Yinit(1:N+1)+Yinit(N+2:2*N+2)+Yinit(2*N+3:3*N+3); %the mass fraction of liquid species in the droplet
    % liquid density constant
    rho_l1 = 995.58-0.1449*Tinit-0.0011*Tinit.^2;%860;
    rho_l2 = 903;
    rho_l3 = 903;
    rho_l = alpha./(Yinit(1:N+1)./rho_l1 + Yinit(N+2:2*N+2)./rho_l2 + Yinit(2*N+3:3*N+3)./rho_l3);
    
    % liquid heat capacity
    cp_l1 = heatcapacityliquid('C8H10',Tinit);
    cp_l2 = heatcapacityliquid('C8H16O2',Tinit);
    cp_l3 = heatcapacityliquid('C8H16O2',Tinit);
    cp_l = (cp_l1.*Yinit(1:N+1) + cp_l2.*Yinit(N+2:2*N+2) + cp_l3.*Yinit(2*N+3:3*N+3))./alpha;
    
    % liquid thermal conductivity
    lambda_l1 = thermalconductivityliquid('C8H10', Tinit);
    lambda_l2 = thermalconductivityliquid('C8H16O2', Tinit);
    lambda_l3 = thermalconductivityliquid('C8H16O2', Tinit);
    lambda_l = alpha./(Yinit(1:N+1)./lambda_l1 + Yinit(N+2:2*N+2)./lambda_l2 + Yinit(2*N+3:3*N+3)./lambda_l3);
    
    % liquid visocisty
    MW1 = 106.16; MW2 = 144.21; MW3 = 405.12;
    Xinit(1:N+1,1)       = Yinit(1:N+1)      /MW1./(Yinit(1:N+1)/MW1 + Yinit(N+2:2*N+2)/MW2 + Yinit(2*N+3:3*N+3)/MW3);
    Xinit(N+2:2*N+2,1)   = Yinit(N+2:2*N+2)  /MW2./(Yinit(1:N+1)/MW1 + Yinit(N+2:2*N+2)/MW2 + Yinit(2*N+3:3*N+3)/MW3);
    Xinit(2*N+3:3*N+3,1) = Yinit(2*N+3:3*N+3)/MW3./(Yinit(1:N+1)/MW1 + Yinit(N+2:2*N+2)/MW2 + Yinit(2*N+3:3*N+3)/MW3);
    viscosity_l1 = viscosityliquid('C8H10',Tinit);
    viscosity_l2 = viscosityliquid('C8H16O2',Tinit);
    viscosity_l3 = viscosityliquid('C8H16O2',Tinit);
    viscosity_l = viscosity_l1.^Xinit(1:N+1) .* viscosity_l2.^Xinit(N+2:2*N+2) .* viscosity_l3.^Xinit(2*N+3:3*N+3);
   
end
