close all;
clear all;

load('gestationdata.mat');
N=100;%number of grid points
t=(1:N)/N;%grid points for scaled response

Y1=DDE_A;
Y2=TRIGLYC;
X=GEST_DEL;
Y=Y1;
ind=find(X<=45);
X=X(ind);%response
Y=Y(ind);%predictor 1
Y1=Y2(ind);%predictor 2
XX=X;

n=length(X);
A=min(X)- (std(X)/sqrt(n));%estimated lower bound for support
B=max(X)+ (std(X)/sqrt(n));%estimated upper bound for support
X=(X-A)/(B-A);%scaling the response

t1=t.*(B-A) + A;%grid for UNscaled data

T = 2*pi*(1:N)/N;


tt0=(1:N)/N;

M=zeros(n,1);
for i=1:n
    g(i)=locallinearregbivariate(XX,Y,Y1,[Y(i,1) Y1(i,1)]);
    err(i)=XX(i)-g(i);
    M(i)=g(i);
    M(i)=(M(i)-A)/(B-A);
end
s1=std(err)/(B-A);
%%% Nonparametric part starts here
%%fourier basis%%
K = 3;
Phi(1:K,:) = sqrt(2)*sin((1:K)'*T);
Phi(K+1:2*K,:) = sqrt(2)*cos((1:K)'*T);
m = size(Phi,1);
%Location of the predictors
y_dde=[quantile(Y,0.1) quantile(Y,.5) quantile(Y,.8) quantile(Y,.95) quantile(Y,.99)];
y_gly=[quantile(Y1,0.1) quantile(Y1,.5) quantile(Y1,.8) quantile(Y1,.95) quantile(Y1,.99)];
y2=y_gly(5); %fixing glyceride location(change according to the problem)
[test11,~,h1]=ksdensity(Y,y_dde);
[test12,~,h2]=ksdensity(Y1,y_gly);
hh=harmmean([h1,h2]);
figure(1);
for j=1:5
    y1=y_dde(j);
    options = optimset('MaxFunEvals',30000,'MaxIter',90000);
    h=hh/sqrt(test11(j)*test12(5));%CHANGE EVERY TIME according to the problem
    strt=zeros(1,m);
    d=fminsearch(@(c)formwteddensityregressionfromMh(c,X,[Y Y1],[y1 y2],N,M,Phi,t,s1,h),strt,options);
        
  d1=d(1:m);
  m1=locallinearregbivariate(XX,Y,Y1,[y1 y2]);
  m1=(m1-A)/(B-A);
  fp=normpdf(t,m1,s1);fp=fp/(sum(fp)/N);
  gamEst = FormGammaFromC(d1,Phi);
  gamDot = gradient(gamEst,1/N);
  fn1 = interp1(t, fp, (t(end)-t(1)).*gamEst + t(1)).*gamDot;
  fn1=fn1/(sum(fn1)/N);%warped density estimate for scaled responses
  fn1=fn1/(B-A);%scaling back the density estimate
  plot(t1,fn1);
  hold on;
  ff(j,:)=fn1;
end