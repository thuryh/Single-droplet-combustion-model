function Tresult = heattransfer(t,T0,rs,alpha_d,lambda_d,reactionrate)
    global N xn dx
    global T_inf Y_O_inf Qc nu cp Bm Qv lambda_g lambda_bd drs2dt
    % the equation is discreted into Ax=BB and solved by LU decomposition.
    % listx, listy indicate the index of the Matrix A where the content is indicated by the list.
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
        BB(i,1) = [T0(i-1),T0(i),T0(i+1)]*[0.25*a(i)*h1+0.5*b(i)*h2, 1-b(i)*h2, -0.25*a(i)*h1+0.5*b(i)*h2]'-reactionrate(i)*(t(end)-t(1))*0;
    end    
    listx = [listx,N+1,N+1];
    listy = [listy,N,  N+1];
    list  = [list, -1, 1-c1*dx];
    BB(N+1,1)= c0*dx;
    
    A = sparse(listx,listy,list,N+1,N+1);
    [L2,U,P] = lu(A);
    Tresult = U\(L2\(P*BB));
end