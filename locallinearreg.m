function g = locallinearreg(Y,Z,x);
n=size(Z,1);
r=size(Z,2);
[a,pts,h]=ksdensity(Z);
z2=normpdf((x-Z)/h,0,1);
W=diag(z2);
e1=[1;zeros(r,1)];
X=cat(2,ones(n,1),x-Z);
g=e1'*inv(X'*W*X)*X'*W*Y;

