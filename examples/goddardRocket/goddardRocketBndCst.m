function output = goddardRocketBndCst(t0,tf,x0,xf,auxdata)
iphase = auxdata.iphase;
if iphase == 1 || iphase == 3
    output = tf - t0;
elseif iphase == 2
    h = x0(1);
    v = x0(2);
    m = x0(3);
    D = auxdata.dragk.*(v.^2).*exp(-h/auxdata.H);
    output = [m*auxdata.g0-(1+v/auxdata.c).*D,tf-t0];
end

end