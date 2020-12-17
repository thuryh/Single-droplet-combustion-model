function Tresult = heattransfer(t,T0,rs,alpha_d,lambda_d)
    global N xn dx
    global T_inf Y_O_inf Qc nu cp Bm Qv lambda_g lambda_bd drs2dt

    h1 = (t(end)-t(1))/dx;
    h2 = (t(end)-t(1))/dx^2;
    a  = -(drs2dt/2/rs^2.*xn + 2*alpha_d/rs^2./xn);
    b  = alpha_d/rs^2;
    c0  = (T_inf+Y_O_inf*Qc/nu/cp)/Bm*log(1+Bm)*lambda_bd/lambda_d(N+1) - Qv/cp*log(1+Bm)*lambda_g/lambda_d(N+1);
    c1  = -log(1+Bm)/Bm*lambda_bd/lambda_d(N+1);
    listx = [];
    listy = [];
    list  = [];
    
    listx = [listx,1,1];
    listy = [listy,1,2];
    list  = [list, 1,-1];
    BB(1,1) = 0;

    for i=2:1:N        
        listx = [listx, i, i, i];
        listy = [listy, i-1, i, i+1];
        list  = [list,  -0.25*a(i)*h1-0.5*b(i)*h2, 1+b(i)*h2, 0.25*a(i)*h1-0.5*b(i)*h2];
        BB(i,1) = [T0(i-1),T0(i),T0(i+1)]*[0.25*a(i)*h1+0.5*b(i)*h2, 1-b(i)*h2, -0.25*a(i)*h1+0.5*b(i)*h2]';
    end    
    listx = [listx,N+1,N+1];
    listy = [listy,N,N+1];
    list  = [list, -1, 1-c1*dx];
    BB(N+1,1)= c0*dx;
    
    A = sparse(listx,listy,list,N+1,N+1);
    [L2,U,P] = lu(A);
    Tresult = U\(L2\(P*BB));
%     C1 = drs2dt/(2*rs^2);
%     C2 = alpha_d/rs^2;
%     %dTdxs  = ((T_inf-Ts+Y_O_inf*Qc/nu/cp)/Bm-Qv/cp)*log(1+Bm)*lambda_g/lambda_l;%((T_inf*cp-Ts*cp+Qc*Y_O_inf/nu)/Bm - Qv)*log(1+Bm)*lambda_g/lambda_l/cp;
%     lambda_l = alpha_d(end)*rho_d(end)*cp_d(end);
%     B1 = ((T_inf+Y_O_inf*Qc/nu/cp)/Bm-Qv/cp)*log(1+Bm)*lambda_g/lambda_l;
%     B0 = 1/Bm*log(1+Bm)*lambda_g/lambda_l;
%     h1 = (t(end)-t(1))/dx;
%     h2 = (t(end)-t(1))/dx^2;
%     listx = [];
%     listy = [];
%     list = [];
%     for i=2:1:N
%         listx = [listx,i,i,i];
%         listy = [listy,i-1,i,i+1];
%         list = [list, 1/2*h1*(C1*xn(i)+2*C2(i)/xn(i))-h2*C2(i), 1+2*C2(i)*h2, -1/2*h1*(C1*xn(i)+2*C2(i)/xn(i))-h2*C2(i)];
%         b(i,1) = rho_d(i)*cp_d(i)*T0(i);
%     end
%     listx = [listx,1,1];
%     listy = [listy,1,2];
%     list = [list,1,-1];
%     b(1,1) = 0;
%     
%     listx = [listx,N+1,N+1,N+1];
%     listy = [listy,N-1,N,N+1];
%     list = [list,1/2/cp_d(end-2)/rho_d(end-2),-2/cp_d(end-1)/rho_d(end-1),(3/2+B0*dx)/cp_d(end)/rho_d(end)];
%     b(N+1) = dx*B1;
%     
%     A = sparse(listx,listy,list,N+1,N+1);
%     [L2,U,P] = lu(A);
%     Tresult = U\(L2\(P*b))./rho_d./cp_d;
end