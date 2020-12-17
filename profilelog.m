
Time_h = [];
rs_h = [];
Tl_h = [];
Tg_h = [];
Y1_h = [];
Y2_h = [];

for i = 0:1:9
    j = 1 - 0.1*i;
    k = min(find(rs_history.^2/rs0^2 <= j));
    Time_h = [Time_h;Time_history(k)];
    rs_h = [rs_h;rs_history(k)];
    Tl_h = [Tl_h;Tresult_history(k,:)];
    Tg_h = [Tg_h;Tgresult_history(k,:)];
    Y1_h = [Y1_h;Y1result_history(k,:)];
    Y2_h = [Y2_h;Y2result_history(k,:)];
    
%     epsilon = YFresult_history(k,1)/YFresult_history(k,3);
%     nu = epsilon*nu1 + (1 - epsilon)*nu2;
%     YFs = YFresult_history(k,3);
%     rf = rf_history(k);
%     Bm = B_history(k,3);
%     
%     YF(find(r< rf)) = -Y_O_inf/nu + (YFs+Y_O_inf/nu)*(exp(-1./r(find(r<rf))*log(1+Bm))-1)/(1/(1+Bm)-1);
%     YF(find(r>=rf)) = zeros(1,length(find(r>=rf)));
%     YO(find(r< rf)) = zeros(1,length(find(r< rf)));
%     YO(find(r>=rf)) = Y_O_inf - (YFs*nu+Y_O_inf)*(exp(-1./r(find(r>=rf))*log(1+Bm))-1)/(1/(1+Bm)-1);
%     YP(find(r< rf)) = (nu+1)/nu*Y_O_inf - (nu+1)/nu*Bm/(1+Bm)*Y_O_inf*(exp(-1./r(find(r< rf))*log(1+Bm))-1)/(1/(1+Bm)-1);
%     YP(find(r>=rf)) = ((nu+1)*YFs+(1+nu)/nu/(1+Bm)*Y_O_inf)*(exp(-1./r(find(r>=rf))*log(1+Bm))-1)/(1/(1+Bm)-1);
%     YN              = Y_N_inf - Bm/(1+Bm)*Y_N_inf*(exp(-1./r*log(1+Bm))-1)/(1/(1+Bm)-1);
%     
%     [MF_CO2,MF_H2O] = GasProducts(epsilon,ST11,ST12,ST21,ST22);
%     YCO2 = YP.*MF_CO2;
%     YH2O = YP.*MF_H2O;
    
    
end
Tl_h = Tl_h';
Tg_h = Tg_h';
Y1_h = Y1_h';
Y2_h = Y2_h';

save('h_10.mat','Time_h','rs_h','Tl_h','Tg_h','Y1_h','Y2_h')