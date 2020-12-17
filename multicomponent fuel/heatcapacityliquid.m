function cpl = heatcapacityliquid(name,T)
    global MW_F1 MW_F2
    if strcmp(name,'C8H16O2')
        a = [3.8452,2.7972,-0.42867,29.24600,-27.587];
        b = [-0.33997,-0.054967,0.93805,3.42610,-0.16485];
        d = [0.19489,0.10679,0.0029498,-2.8920,3.05310];
        n = [2,4,1,1,1];
        A = sum(a.*n);
        B = sum(b.*n);
        D = sum(d.*n);
        MW = MW_F2;
    end
    if strcmp(name,'C8H10')
        a = [3.8452,1.9570,3.6968];
        b = [-0.33997,-0.31938,-1.6037];
        d = [0.19489,0.11911,0.55022];
        n = [2,2,4];
        A = sum(a.*n);
        B = sum(b.*n);
        D = sum(d.*n);
        MW = MW_F1;
    end
    cpl = 8.314*(A+B*T/100+D*(T/100).^2)/MW;
end