close all;
figure(1)
imagesc(xn,Time_history/rs0^2/1e6/4,Tresult_history,[300,T_boil_F2]);colormap('jet');colorbar;
hold on;contour(xn,Time_history/rs0^2/1e6/4, Tresult_history, [300,T_boil_F1,T_boil_F2],'k-');
%hold on;contour(xn,Time_history/rs0^2/1e6/4, Tresult_history./Tb', [0.5,1],'b-');
xlabel('Nondimensional radial position, R=r/r_s');ylabel('Droplet lifetime, t/d_{s0}^2 \mus/\mum^2');
set(gca, 'FontSize',14,'FontName','Arial');set(get(colorbar,'title'),'string','T (K)');
figure(2)
imagesc(xn,Time_history/rs0^2/1e6/4,Y1result_history,[0,1]);colormap('jet');colorbar;
xlabel('Nondimensional radial position, R=r/r_s');ylabel('Droplet lifetime, t/d_{s0}^2 \mus/\mum^2');
set(gca, 'FontSize',14,'FontName','Arial');set(get(colorbar,'title'),'string','Y');

%load m-xylene-2-EHA-Le=10-pureoxygen-case.mat
%global Qv1 Qv2 MW_F1 MW_F2 Rg
%Tb = [];
% for i=1:1:length(Time_history)
%     Tb = [Tb,Tboilpoint([Y1result_history(i,:)';Y2result_history(i,:)'],T_boil_F1,T_boil_F2)];
%     i
% end

output1 = [xn,Y1result_history(400,:)',Y2result_history(400,:)',Y1result_history(1200,:)',Y2result_history(1200,:)',Y1result_history(2400,:)',Y2result_history(2400,:)'];
