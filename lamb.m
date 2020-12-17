function lambda = lamb(Y,X,T)
global Patm Rg kb Na
N = 6;
% C8H10 C8H16O2 CO2 O2 N2 H2O

MW = [106;144;44;32;28;18]/1e3;
sigma_k = [6.51;7.55;3.76;3.46;3.62;2.60].*1e-10;
C_v =[specheat('C8H10',T);specheat('C8H16O2',T);specheat('CO2',T);specheat('O2',T);specheat('N2',T);specheat('H2O',T)]- Rg;% refer to specific heat; C_v = C_p - R
Z_rot_0 = [0.00;0.00;2.10;3.80;4.00;4.00];
epsilon = [4003.83;4586.38;2028.74;892.98;810.91;4759.22]/Na;
myu_k = [0.00;0.00;0.00;0.00;0.00;1.84]; % A s m

T_s = kb*T./epsilon;
Omega_D = 1.06036./T_s.^0.15610 + 0.19300./exp(0.47635*T_s) + 1.03587./exp(1.52996*T_s) + 1.76474./exp(3.89411*T_s);
rho = Patm*MW/Rg/T;
D_kk = 3/8*sqrt(pi*kb^3*T^3./MW*Na)/Patm/pi./sigma_k.^2./Omega_D;
    
delta_k = 1/2*myu_k.^2/epsilon./sigma_k.^3;
Omega_22 = 1.16145./T_s.^0.14874 + 0.52487./exp(0.7732*T_s) + 2.16178./exp(2.43787*T_s)-6.435e-4*T_s.^0.14874.*sin(18.0323*T_s.^(-0.7683)-7.27371);
        
eta_k = 5/16*sqrt(pi*MW/Na*kb*T)/pi./sigma_k.^2./Omega_22;

C_v_trans = 1.5*Rg;
C_v_rot = [1.5, 1.5, 1, 1, 1, 1.5]'*Rg;
C_v_vib = C_v - [3, 3, 2.5, 2.5, 2.5, 3]'*Rg;
Z_rot = Z_rot_0.*F(298,epsilon)./F(T,epsilon);
    
A = 5/2 - rho.*D_kk./eta_k;
B = Z_rot + 2/pi*(5/3*C_v_rot/Rg + rho.*D_kk./eta_k);
    
f_trans = 5/2*(1 - 2/pi*C_v_rot./C_v_trans.*A./B);
f_rot = rho.*D_kk./eta_k.*(1 + 2/pi*A./B);
f_vib = rho.*D_kk./eta_k;
    
lambda_k = eta_k./MW.*(f_trans.*C_v_trans + f_rot.*C_v_rot + f_vib.*C_v_vib);
lambda = 1/2*(sum(X'.*lambda_k) + 1./sum(X'./lambda_k));
end

function F = F(T,epsilon)
global kb
F = 1 + pi^(3/2)/2*(epsilon/kb/T).^(1/2) + (pi^2/4 + 2)*(epsilon/kb/T) + pi^(3/2)*(epsilon/kb/T).^(3/2);
end

