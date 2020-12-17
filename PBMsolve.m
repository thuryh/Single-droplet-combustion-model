function dPBM = PBMsolve(t,PBM,Sp,rs,v,rho_g,T,drs2dt)
    global alpha MW_g kb Na rho_p vp0
    global r dr
    
    M = length(r);
    psi = PBM(1:M);
    n   = PBM(M+1:end);
    dpsi= zeros(M,1);
    dn  = zeros(M,1);
    [vp,dp] = calvpdp(psi,n);
    
    D = 3./(2*(1+pi*alpha/8)*dp.^2.*rho_g').*sqrt(MW_g'*kb.*T'/(2*pi)/Na); % particle diffusivity, m2 s-1
    beta = (3/(4*pi))^(1/6) * (6*kb*T'/rho_p).^0.5 * sqrt(2)^2.5 .* vp.^(1/6);%particle coagulation rate, m3 s-1
    
%% psi equations
% fatal convection of psi
    rc = (r(1:end-1) + r(2:end))/2; % 1.5,2.5,...rN-0.5
    mconvpsi = diff(psi)/dr.*rc'; %1.5,2.5...rN-0.5 * 2,3,4...rN
% flux
    Dc = (D(1:end-1)+D(2:end))/2;
    fluxpsi1 = (max(v(1:end-1)',0).*psi(1:end-1) + min(v(2:end)',0).*psi(2:end)) * rs;
    fluxpsi2 = -Dc.*diff(psi)/dr;
    fluxpsi = (fluxpsi1+fluxpsi2).*rc'.^2;
    dfluxpsi = [diff([fluxpsi;0])/dr./r(2:end)'.^2];
% final equation
    dpsi(2:end) = Sp(2:end)' + mconvpsi*drs2dt/2/rs^2 - dfluxpsi/rs^2 ;
    dpsi(1) = 0;
    dpsi(end) = 0;
%% n equation
% fatal convection of n
    mconvn = diff(n)/dr.*rc'; %1.5,2.5...rN-0.5 * 2,3,4...rN
% flux
    fluxn1 = (max(v(1:end-1)',0).*n(1:end-1) + min(v(2:end)',0).*n(2:end)) * rs;
    fluxn2 = -Dc.*diff(n)/dr;
    fluxn = (fluxn1+fluxn2).*rc'.^2;
    dfluxn = diff([fluxn;0])/dr./r(2:end)'.^2;
%final equation
    dn(2:end) = Sp(2:end)'./vp0 - 1/2*beta(2:end).*n(2:end).^2 - mconvn*drs2dt/2/rs^2 - dfluxn/rs^2 ;
    dn(1) = 0;
    dn(end) = 0;

%% output
    dPBM = [dpsi;dn];
end