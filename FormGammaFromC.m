function gam = FormGammaFromC(c,p)

N = size(p,2);
v = c*p;

nc = norm(c,'fro');

if (nc==0)
    q=ones(1,N);
else
q = cos(nc)*ones(1,N)+sin(nc)*v/nc;
end

gam = cumsum(q.*q)/N;
gam = (gam-gam(1))/(gam(end)-gam(1));

