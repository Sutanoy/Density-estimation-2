close all;
clear all;

load('gestationdata.mat');
N=100;
t=(1:N)/N;

Y1=DDE_A; %DDE dose
Y2=TRIGLYC;%triglyceride level
X=GEST_DEL; %Gestation age at delivery
Y=Y1; %predictor. Change according to the requirements of the analysis.
ind=find(X<=45);
X=X(ind); %response
Y=Y(ind); %predictor
XX=X;

n=length(X);
A=min(X)- (std(X)/sqrt(n));%lower support boundary estimate
B=max(X)+ (std(X)/sqrt(n));%upper boundary estimate
X=(X-A)/(B-A); %scaling
t1=t.*(B-A) + A;

T = 2*pi*(1:N)/N;


tt0=(1:N)/N;

M=zeros(n,1);
for i=1:n
    g(i)=locallinearreg(XX,Y,Y(i,1));
    err(i)=XX(i)-g(i);
    M(i)=g(i);
    M(i)=(M(i)-A)/(B-A);
end
s1=std(err)/(B-A);
%%% Nonparametric part starts here
%%fourier basis%%
K = 3; %As with simulated examples, 6 basis elements used
Phi(1:K,:) = sqrt(2)*sin((1:K)'*T);
Phi(K+1:2*K,:) = sqrt(2)*cos((1:K)'*T);
m = size(Phi,1);

y=[quantile(Y,0.1) quantile(Y,.5) quantile(Y,.8) quantile(Y,.95) quantile(Y,.99)];
%^The locations at which we compute the conditional density estimation
[test1,test2,h]=ksdensity(Y,y);%The naive density estimate at the locaions and bandwidth
figure(1);clf;
for j=1:length(y)
   options = optimset('MaxFunEvals',30000,'MaxIter',90000);
        h1=h/sqrt(test1(j));
        strt=zeros(1,m);
        d=fminsearch(@(c)formwteddensityregressionfromMh(c,X,Y,y(j),N,M,Phi,t,s1,h1),strt,options);
        d1=d(1:m);
  m1=locallinearreg(XX,Y,y(j));
  m1=(m1-A)/(B-A);
  fp=normpdf(t,m1,s1);fp=fp/(sum(fp)/N);
  gamEst = FormGammaFromC(d1,Phi);
  gamDot = gradient(gamEst,1/(N));
  fn1 = interp1(t, fp, (t(end)-t(1)).*gamEst + t(1)).*gamDot;
  fn1=fn1/(sum(fn1)/N);
  fn1=fn1/(B-A);
  plot(t1,fn1);
  hold on;
  ff(j,:)=fn1;
end