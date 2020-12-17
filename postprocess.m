% Tcenter = Tresult_history(:,1);
Tsurface= Tresult_history(:,end);

% Y1center = Y1result_history(:,1);
% Y1surface= Y1result_history(:,end);
% 
% Y2center = Y2result_history(:,1);
% Y2surface= Y2result_history(:,end);
% 
% rs2_history = rs_history.^2;
% 
% M = 1 - rs2_history/rs0^2;
% 
%     %% For saving dp, psi, and workspace data
% 
%     figure(1);imagesc(r(2:end),Time_history(1:end-1)*1e6,psiresult_history(1:end-1,2:end));
%     xlabel('Nondimensional Radius, r/r_s');ylabel('Time, t (\mu s)');
%     set(gca,'YDir','reverse','FontSize',14,'FontName','Arial');
%     colorbar;
%     saveas(gcf,['psi_r=',num2str(rs0*1e6),',Le=',num2str(Le_l1),',eha=',num2str(Y0(N+2)*1e2),'.bmp']);
%     
%     
%     
    figure(2);imagesc(r(2:end),Time_history(1:end-1)*1e6,dpresult_history(1:end-1,2:end)*1e9);
    xlabel('Nondimensional Radius, r/r_s');ylabel('Time, t (\mu s)');
    set(gca,'YDir','reverse','FontSize',14,'FontName','Arial');
    colorbar;
    saveas(gcf,['dp_r=',num2str(rs0*1e6),',Le=',num2str(Le),',Al(NO3)3=',num2str(Y0(N+2)*1e2),'.bmp']);
%     
    save(['r=',num2str(rs0*1e6),',Le=',num2str(Le),',Al(NO3)3=',num2str(Y0(N+2)*1e2),',Tinf =',num2str(T_inf),'.mat']);
%     
    %% For validation of Fuel mass loss
    
    Sigmadm = zeros(length(dmdt_history),1);
    Sigmadm(1) = dmdt_history(1)*Time_history(1);
    for i = 2:length(Sigmadm)
        Sigmadm(i) = dmdt_history(i).*(Time_history(i)-Time_history(i-1)) + Sigmadm(i-1);
    end
    EtOHleftb = 4/3*pi*rs0^3*rho_l*Y0(1) - Sigmadm;
    
    massleft = 4/3*pi*rs0^3*rho_l - Sigmadm;
    
    EtOHleftc = zeros(length(Time_history)-1,1);
    for i = 1:length(EtOHleftc)
        EtOHleftc(i) = sum(4*pi*rs_history(i)^3*((xn(1:end-1)+xn(2:end))/2).^2.*((Y1result_history(i,1:end-1)+Y1result_history(i,2:end))/2)'*rho_l*dx);
    end
    plot(Time_history, EtOHleftb,Time_history(1:end-1), EtOHleftc);
    
    
    %% For comparison of Temperature - Time
    Tsurface= Tresult_history(:,end);
    % Al(NO3)3 decomposition
    A = 7.58E15;
    E = 1.55E5;
    t_chemi = 1./(A*exp(-E./Tsurface/Rg));
    
    [tpoint,Point] = max(min(t_chemi(1:end-1),Time_history));
    
%     Tpoint = E/Rg/log(A*tpoint);
    Tpoint = Tsurface(Point);
    Rpoint = rs_history(Point);
    TimeRes = massleft(Point)/dmdt_history(Point);
    
    fprintf('t = %1.4e s, T = %1.4e K, rs = %1.4e m, Res Time = %1.4e s\n',tpoint,Tpoint,Rpoint,TimeRes);
    
    
    Nondimtime = Time_history/Time_history(end);
    Nondimrs2 = rs_history.^2/rs0^2;
    [Dpm,Posi] = max(dpresult_history(:,5));
    fprintf('rs = %1.4e m, Dpm = %1.4e m, Psi = %1.4e, Ypre =  %1.4e  \n',rs_history(Posi),Dpm,psiresult_history(Posi,5), YFresult_history(Posi,2));
    clear;
    Tsurface= Tresult_history(:,end);
    Point = length(Tsurface(find(Tsurface<=T_boil_F1)));
%     Tsurface(Point);
    Time_history(Point)
    TimeRes = massleft(Point)/dmdt_history(Point)
    
    imagesc(r,Time_history,Y2result_history(1:end-1,:));
    xlabel('Nondimensional Radius, r/r_s');ylabel('Time, t (\mu s)');
    set(gca,'YDir','reverse','FontSize',14,'FontName','Arial');
    colorbar;
    
    Tsurface= Tresult_history(:,end);
    Y2surface= Y2result_history(:,end);
plot(Y2result_history(:,end),Tsurface)


    Tsurface= Tresult_history(:,end);
    Point = length(Tsurface(find(Tsurface<=T_boil_F1)));
    Point = length(Tsurface);
color = [];
for i = 1:20
    color = [color;i/20,0.2,0.4];
end
figure;
for i=1:20
    j = floor(Point/20)*i;
    plot(xn*rs_history(j),Y2result_history(j,:),'color',color(i,:));
%     leg_str{i-1}=[num2str(SS(i-1)),'mg/L'];
    hold on
end
xn

plot(Time_history(1:end-5),dmdt_history(1:end-5));ylabel('Mass Loss Rate, kg/s');xlabel('Time, t (\mu s)');

%%%%%%%%%%%%%%%20200217
    [rsize,~] = size(Time_history);

    figure
    imagesc(xn,Time_history(1:rsize-5)*1e6,Y1result_history(1:rsize-5,:));
    xlabel('Nondimensional Radius, r/r_s');ylabel('Time, t (\mu s)');
    set(gca,'YDir','reverse','FontSize',14,'FontName','Arial');
    colorbar;
    title('Y1');
    saveas(gcf,['Y1','.bmp'])
    hold on
    
    figure
    imagesc(xn,Time_history(1:rsize-5)*1e6,Y2result_history(1:rsize-5,:));
    xlabel('Nondimensional Radius, r/r_s');ylabel('Time, t (\mu s)');
    set(gca,'YDir','reverse','FontSize',14,'FontName','Arial');
    colorbar;
    title('Y2');
    saveas(gcf,['Y2','.bmp'])
    hold on
    
    figure
    imagesc(xn,Time_history(1:rsize-5)*1e6,1 - phiresult_history(1:rsize-5,:));
    xlabel('Nondimensional Radius, r/r_s');ylabel('Time, t (\mu s)');
    set(gca,'YDir','reverse','FontSize',14,'FontName','Arial');
    title('\phi');
    colorbar;
    saveas(gcf,['phi','.bmp'])
    hold on
    
    figure
    plot(Time_history(1:rsize-5)*1e6,(rs_history(1:rsize-5)*1e6).^2);
    ylabel('Radius^2, r^2 (\mu m^2)');xlabel('Time, t (\mu s)');
    title('r_2');
    saveas(gcf,['rs2','.bmp'])
    
    figure
    plot(Time_history(1:rsize-5)*1e6,(Tresult_history(1:rsize-5,end)));
    ylabel('Temperature, T (K)');xlabel('Time, t (\mu s)');
    title('Liquid Phase Temp');
    saveas(gcf,['Ts','.bmp'])
    hold on
    
    figure
    imagesc(r,Time_history(1:rsize-5)*1e6,Tgresult_history(1:rsize-5,:));
    xlabel('Nondimensional Radius, r/r_s');ylabel('Time, t (\mu s)');
    set(gca,'YDir','reverse','FontSize',14,'FontName','Arial');
    colorbar;
    title('Gas Phase Temp');
    saveas(gcf,['Tg','.bmp'])
    hold on
    
    figure
    imagesc(xn,Time_history(1:rsize-5)*1e6,Tresult_history(1:rsize-5,:));
    xlabel('Nondimensional Radius, r/r_s');ylabel('Time, t (\mu s)');
    set(gca,'YDir','reverse','FontSize',14,'FontName','Arial');
    colorbar;
    title('Absolute Temp');
    saveas(gcf,['T_abs','.bmp'])
    hold on
    
    figure
    imagesc(xn,Time_history(1:rsize-5)*1e6,Tresult_history(1:rsize-5,:)-T_boil_F1);
    xlabel('Nondimensional Radius, r/r_s');ylabel('Time, t (\mu s)');
    set(gca,'YDir','reverse','FontSize',14,'FontName','Arial');
    colorbar;
    title('Companion Temp');
    saveas(gcf,['T_com','.bmp'])
    hold on
%     [rsize,~] = size(Time_history);
    figure
    imagesc(xn,Time_history(1:rsize-5)*1e6,Y2result_history(1:rsize-5,:)/0.2);
    xlabel('Nondimensional Radius, r/r_s');ylabel('Time, t (\mu s)');
    set(gca,'YDir','reverse','FontSize',14,'FontName','Arial');
    colorbar;
    title('S');
    saveas(gcf,['S','.bmp'])
    hold on
    
    figure
    imagesc(xn,Time_history(1:rsize-5)*1e6,ppsiresult_history(1:rsize-5,:));
    xlabel('Nondimensional Radius, r/r_s');ylabel('Time, t (\mu s)');
    set(gca,'YDir','reverse','FontSize',14,'FontName','Arial');
    colorbar;
    title('\Psi');
    saveas(gcf,['Psi','.bmp'])
    hold on
    
    figure
    plot(Time_history(1:rsize-5)*1e6,(rs_history(1:rsize-5).*rf_history(1:rsize-5)*1e6));
    ylabel('Radius, r (\mu m)');xlabel('Time, t (\mu s)');
    title('rf');
    
    figure
    plot(Time_history(1:rsize-5)*1e6,(rs_history(1:rsize-5)*1e6));
    ylabel('Radius, r (\mu m)');xlabel('Time, t (\mu s)');
    title('r');
    
    figure
    plot(Time_history(1:rsize-5)*1e6,rf_history(1:rsize-5));
    ylabel('Radius, r');xlabel('Time, t (\mu s)');
    title('rf');
    
    figure
    plot(Time_history(1:rsize-5)*1e6,(rs_history(1:rsize-5)/rs0).^2);
    ylabel('Radius, r (\mu m)');xlabel('Time, t (\mu s)');
    title('r');
    
    figure
    semilogx(Time_history(1:rsize-5)*1e6,Tresult_history(1:rsize-5,end));
    ylabel('Temperature, K');xlabel('Time, t (\mu s)');
    title('Surface Temp');
    
    
    figure
    plot(Time_history(1:rsize-5)*1e6,Y1result_history(1:rsize-5,end));
    ylabel('M - Xylene Concentration, /');xlabel('Time, t (\mu s)');
    title('Surface Y');