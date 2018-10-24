function output = goddardRocketBndObj(t0,tf,x0,xf,auxdata)

iphase = auxdata.iphase;
if iphase == 1 || iphase == 2
    output = zeros(size(tf));
else 
    output = -xf(1);
end

end