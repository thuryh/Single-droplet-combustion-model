function phiresult = particletransfer(t,phiinit,rho,rs,DS,omega)

    global N xn dx
    global lambda_g cp drs2dt
    global Bm rho_s MW_S MW_Pr MW_F2
    % the equation is discreted into Ax=BB and solved by LU decomposition.
    % listx, listy indicate the index of the Matrix A where the content is indicated by the list.
    dt = t(end)-t(1);
    h1 = dt/dx;
    h2 = dt/dx^2;
    
    %% solving phi
    a  = -(drs2dt/2/rs^2.*xn + 2*DS/rs^2./xn);
    b  = DS/rs^2;
    B0 = log(1+Bm)*lambda_g/cp/DS(N+1)/rho(end);

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
        BB(i,1) = [phiinit(i-1),phiinit(i),phiinit(i+1)]*[0.25*a(i)*h1+0.5*b(i)*h2, 1-b(i)*h2, -0.25*a(i)*h1+0.5*b(i)*h2]'+omega(i)*dt/rho_s*(MW_Pr-2*MW_F2)/MW_Pr;
    end
    listx = [listx,N+1,N+1,N+1];
    listy = [listy,N-1,N,  N+1];
    list  = [list, 1/2,-2, 3/2-B0*dx];
    BB(N+1) = 0;
    
    A = sparse(listx,listy,list,N+1,N+1);
    [L2,U,P] = lu(A);
    phiresult = U\(L2\(P*BB));   
end