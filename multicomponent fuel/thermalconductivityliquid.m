function lambda = thermalconductivityliquid(name, T)
global T_boil_F1 T_boil_F2 MW_F1 MW_F2
    if strcmp(name,'C8H16O2')
        Tb = 214 + 273.15;
        %CH3 CH2 CH COOH
        Nktck = [0.0141, 0.0189, 0.0164, 0.0791];
        item = [2,4,1,1];
        Tc = Tb*(0.584+0.965*sum(Nktck.*item)-(sum(Nktck.*item))^2)^(-1);
        Astar = 0.00319;
        alpha = 1.2;
        beta = 0.5;
        gamma = 0.167;
        A = Astar*Tb^alpha/(MW_F2*1e3)^beta/Tc^gamma;
        Tr = T/Tc;
    end
    if strcmp(name,'C8H10')
        Tb = 139 + 273.15;
        %CH3 CH= C=
        Nktck = [0.0141, 0.0082, 0.0143];
        item = [2,4,1];
        Tc = Tb*(0.584+0.965*sum(Nktck.*item)-(sum(Nktck.*item))^2)^(-1);
        Astar = 0.00319;
        alpha = 1.2;
        beta = 0.5;
        gamma = 0.167;
        A = Astar*Tb^alpha/(MW_F1*1e3)^beta/Tc^gamma;
        Tr = T/Tc;
    end
    lambda = A*(1-Tr).^0.38./(Tr).^(1/6);

end