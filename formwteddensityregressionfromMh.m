function L=formwteddensityregressionfromMh(c,X,Y,x0,N,M,Phi,t,s1,h)
z=c;

gam0 = FormGammaFromC(z,Phi);
gam = (gam0-gam0(1))/(gam0(end)-gam0(1));
gamDot = gradient(gam,1/N);
gamatx=interp1(t,gam,X,'linear','extrap');
gamDotx=interp1(t,gamDot,X,'linear','extrap');
fpnorm=normcdf(1,M,s1)-normcdf(0,M,s1);
fpnew=normpdf(gamatx,M,s1)./fpnorm;
diff=Y-ones(length(X),1)*x0;
normdiff=sqrt(sum(diff.^2,2));
wt=normpdf(normdiff/h,0,1);
wt=wt/sum(wt);
%L=-sum(log(fpnew.*gamDotx.*wt));
ind=find(wt>0.000001);
% ind=find(wt>0);
L=-sum(log(fpnew(ind).*gamDotx(ind).*wt(ind)));