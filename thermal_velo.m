Tch = 1/6*T(1) + 2/3*T(rfp) + 1/6*T(end);
    
YFch = 1/6*YF(1) + 2/3*YF(rfp) + 1/6*YF(end);
YOch = 1/6*YO(1) + 2/3*YO(rfp) + 1/6*YO(end);
YPch = 1/6*YP(1) + 2/3*YP(rfp) + 1/6*YP(end);
YNch = 1/6*YN(1) + 2/3*YN(rfp) + 1/6*YN(end);

dTdr(find(r< rf))  = (Ts - T_inf - Y_O_inf*Qc/nu/cp)*exp((-1./r(find(r< rf))*log(1+Bm)))*log(1+Bm)/(1/(1+Bm)-1)./r(find(r< rf)).^2/rs;
dTdr(find(r>= rf)) = (Ts - T_inf + YFs*Qc/cp)*exp((-1./r(find(r>= rf))*log(1+Bm)))*log(1+Bm)/(1/(1+Bm)-1)./r(find(r>= rf)).^2/rs;
muF1 = 0.01126;% Dortmund Databank
muF2 = 9.7e-6;% 2-Ethylhexanol, MSDS
%Sutherland's formula of viscosity
muO2 = (20 + 273.15 + 125)/(Tch + 125)*(Tch/(273.15 + 20))^(3/2)*20.4*1e-6;
muCO2 = (20 + 273.15 + 240)/(Tch + 240)*(Tch/(273.15 + 20))^(3/2)*14.7*1e-6;
muH2O = (100 + 273.15 + 650)/(Tch + 650)*(Tch/(273.15 + 100))^(3/2)*12.1*1e-6;
muN2 = (20 + 273.15 + 104)/(Tch + 104)*(Tch/(273.15 + 20))^(3/2)*17.6*1e-6;
% Grunberg–Nissan relation
LNmu = log(muF1)*(YFch*epsilon) + log(muF2)*(YFch*(1-epsilon)) + log(muO2)*YOch + (log(muCO2) + log(muH2O))/2*YPch + log(muN2)*YNch;
mu = exp(LNmu);
%nyu = mu./rho_g;
vth = -3.*mu.*dTdr./(4*(1+pi*alpha/8).*T);


