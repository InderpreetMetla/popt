function output = goddardRocketPathCst(t,x,u,auxdata)
iphase = auxdata.iphase;
if iphase == 1 || iphase == 3
    output = [];
elseif iphase == 2
    T = u;
    h = x(:,1);
    v = x(:,2);
    m = x(:,3);
    D = auxdata.dragk.*(v.^2).*exp(-h/auxdata.H);
    voverc = v/auxdata.c; 
    term1 = (auxdata.c^2).*(ones(size(t))+voverc)./...
        (auxdata.H*auxdata.g0)-ones(size(t))-2./voverc;
    term2 = m*auxdata.g0./(ones(size(t))+4./voverc+2./(voverc.^2));
    output = T-D-m*auxdata.g0-term1.*term2;
end
end