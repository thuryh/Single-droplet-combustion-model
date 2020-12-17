function Yresult = masstransfer(t,Y0,phi0,rho,rs,D,omega)
    global N xn dx
    global lambda_g lambda_bd cp cp_bd drs2dt
    global Bm epsilon1 epsilon2 epsilon3 MW_F2 MW_Pr rho_s rho_l
    % the equation is discreted into Ax=BB and solved by LU decomposition.
    % listx, listy indicate the index of the Matrix A where the content is indicated by the list.
    Y10 = Y0(    1:  N+1)./(1-phi0*rho_s./rho);% the mass fraction of species in the liquid phase
    Y20 = Y0(  N+2:2*N+2)./(1-phi0*rho_s./rho);% the mass fraction of species in the liquid phase
    Y30 = Y0(2*N+3:3*N+3)./(1-phi0*rho_s./rho);% the mass fraction of species in the liquid phase
    
    h1 = (t(end)-t(1))/dx;
    h2 = (t(end)-t(1))/dx^2;
    
    %% solving Y10
    a  = -(drs2dt/2/rs^2.*xn + 2*D/rs^2./xn);
    b  = D/rs^2;
    B1 = -log(1+Bm)*epsilon1*lambda_g/cp/rho(end)/D(N+1);
    B0 = -log(1+Bm)*lambda_g/cp/D(N+1)/rho(end);

    listx = [];
    listy = [];
    list  = [];
    
    listx = [listx,1,1];
    listy = [listy,1,2];
    list  = [list, 1,-1];
    BB(1,1) = 0;

    for i=2:1:N
        listx = [listx,i,  i,i  ];
        listy = [listy,i-1,i,i+1];
        list  = [list,  -0.25*a(i)*h1-0.5*b(i)*h2, 1+b(i)*h2, 0.25*a(i)*h1-0.5*b(i)*h2];
        BB(i,1) = [Y10(i-1),Y10(i),Y10(i+1)]*[0.25*a(i)*h1+0.5*b(i)*h2, 1-b(i)*h2, -0.25*a(i)*h1+0.5*b(i)*h2]';
    end
    listx = [listx,N+1,N+1,N+1];
    listy = [listy,N-1,N,  N+1];
    list  = [list, 1/2,-2, 3/2+B0*dx];
    BB(N+1) = dx*B1;
    
    A = sparse(listx,listy,list,N+1,N+1);
    [L2,U,P] = lu(A);
    Y1result = U\(L2\(P*BB));
    
    %% solving Y20
    a  = -(drs2dt/2/rs^2.*xn + 2*D/rs^2./xn);
    b  = D/rs^2;
    B1 = -log(1+Bm)*epsilon2*lambda_g/cp/rho(end)/D(N+1);
    B0 = -log(1+Bm)*lambda_g/cp/D(N+1)/rho(end);

    listx = [];
    listy = [];
    list  = [];
    
    listx = [listx,1,1];
    listy = [listy,1,2];
    list  = [list, 1,-1];
    BB(1,1) = 0;

    for i=2:1:N
        listx = [listx,i,  i,i  ];
        listy = [listy,i-1,i,i+1];
        list  = [list,  -0.25*a(i)*h1-0.5*b(i)*h2, 1+b(i)*h2, 0.25*a(i)*h1-0.5*b(i)*h2];
        BB(i,1) = [Y20(i-1),Y20(i),Y20(i+1)]*[0.25*a(i)*h1+0.5*b(i)*h2, 1-b(i)*h2, -0.25*a(i)*h1+0.5*b(i)*h2]'+2*MW_F2/MW_Pr*omega(i)/rho(i)*(t(end)-t(1));
    end
    listx = [listx,N+1,N+1,N+1];
    listy = [listy,N-1,N,  N+1];
    list  = [list, 1/2,-2, 3/2+B0*dx];
    BB(N+1) = dx*B1;
    
    A = sparse(listx,listy,list,N+1,N+1);
    [L2,U,P] = lu(A);
    Y2result = U\(L2\(P*BB));    
        
    %% solving Y30
    a  = -(drs2dt/2/rs^2.*xn + 2*D/rs^2./xn);
    b  = D/rs^2;
    B1 = -log(1+Bm)*epsilon3*lambda_g/cp/rho(end)/D(N+1);
    B0 = -log(1+Bm)*lambda_g/cp/D(N+1)/rho(end);

    listx = [];
    listy = [];
    list  = [];
    
    listx = [listx,1,1];
    listy = [listy,1,2];
    list  = [list, 1,-1];
    BB(1,1) = 0;

    for i=2:1:N
        listx = [listx,i,  i,i  ];
        listy = [listy,i-1,i,i+1];
        list  = [list,  -0.25*a(i)*h1-0.5*b(i)*h2, 1+b(i)*h2, 0.25*a(i)*h1-0.5*b(i)*h2];
        BB(i,1) = [Y30(i-1),Y30(i),Y30(i+1)]*[0.25*a(i)*h1+0.5*b(i)*h2, 1-b(i)*h2, -0.25*a(i)*h1+0.5*b(i)*h2]'-omega(i)/rho(i)*(t(end)-t(1));
    end
    listx = [listx,N+1,N+1,N+1];
    listy = [listy,N-1,N,  N+1];
    list  = [list, 1/2,-2, 3/2+B0*dx];
    BB(N+1) = dx*B1;
    
    A = sparse(listx,listy,list,N+1,N+1);
    [L2,U,P] = lu(A);
    Y3result = U\(L2\(P*BB));
    %% transforming to the mass fraction of species in the droplet
    Yresult = [(1-Y2result-Y3result).*(1-phi0*rho_s./rho); Y2result.*(1-phi0*rho_s./rho); Y3result.*(1-phi0*rho_s./rho)];
end