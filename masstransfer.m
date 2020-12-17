function Yresult = masstransfer(t,Y0,rho,rs,D)
    global N xn dx
    global lambda_g lambda_bd cp cp_bd drs2dt
    global Bm epsilon
    
    Y10 = Y0(1:N+1);
    Y20 = 1 - Y10;
    
    h1 = (t(end)-t(1))/dx;
    h2 = (t(end)-t(1))/dx^2;
    
    a  = -(drs2dt/2/rs^2.*xn + 2*D/rs^2./xn);
    b  = D/rs^2;
    B1 = -log(1+Bm)*epsilon*lambda_g/cp/rho(end)/D(N+1);
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
    Yresult = [Y1result;1-Y1result];
end