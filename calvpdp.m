function [vp,dp] = calvpdp(psi,n)
global vp0 dp0
    items = find(n.*psi > 0);
    vp = ones(size(psi))*vp0;
    if(~isempty(items))
        vp(items) = psi(items)./n(items);
    end
    dp = (vp*6/pi).^(1/3);
    n(psi==0) = 0;
end
