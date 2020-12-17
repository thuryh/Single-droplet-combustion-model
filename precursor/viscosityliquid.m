function viscosity = viscosityliquid(name, Tinit)
    
    rho_l1 = 995.58-0.1449*Tinit-0.0011*Tinit.^2;
    rho_l2 = 903;
    rho_l3 = 1251;
    A = [];
    if strcmp(name,'C8H10')
        A = -(6.95+0.21*8)+0+0.05;
        B = 275+99*8+20-34;
        rhol = 860;
        MW = 106.16;
    end
    if strcmp(name,'C8H16O2')
        A = -6.95-0.21*8 - 0.15 -0.9;
        B = 275+99*8 + 35 + 770;
        rhol = 903;
        MW = 144.214;
    end
    if strcmp(name,'C16H30O4Sn')
        A = -6.95-0.21*16 - 0.15*2 -0.9*2;
        B = 275+99*16 + 35*2 + 770*2;
        rhol = 1.251e3;
        MW = 405.12;
    end
    viscosity = rhol/1e3*MW.*exp(A+B./Tinit)*0.001;
    
end