clear all; 
close all;

N=100;
t=(1:N)/N;
n=100;

T = 2*pi*(1:N)/N;
clear A;

K = 3;
Phi(1:K,:) = sqrt(2)*sin((1:K)'*T);
Phi(K+1:2*K,:) = sqrt(2)*cos((1:K)'*T);
m = size(Phi,1);
t0=0:0.1:1;
t00=0:0.1:1;
err=zeros(n,1);
%an illustrative simulation example%%%
   Y=rand(n,1);%predictor 1
   Y1=rand(n,1);%predictor 2
   Y11=exp(-Y1);
   u=rand(n,1);

   X=((2*Y -1).^2+laprnd(n,1,0,1)).*(u<Y11) + (u>=Y11).*normrnd((3+2*Y),0.5);
   %response example
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   XXX=X;
   A=min(X)- (std(X)/sqrt(n));%lower boundary of support
   B=max(X)+ (std(X)/sqrt(n));%upper boundary
   X=(X-A)/(B-A); % scaled
   t1=t.*(B-A) + A; %transformed grid points to evaluate final density after rescaling
   
   M=zeros(n,1);
   for i=1:n
    M(i)=locallinearregbivariate(XXX,Y,Y1,[Y(i,1) Y1(i,1)]);
    err(i)=XX(i)-M(i);
    M(i)=(M(i)-A)/(B-A);
   end
s1=std(err)/(B-A);
  [test11,~,h1]=ksdensity(Y,t0);
   [test12,~,h2]=ksdensity(Y1,t00);
   hh=harmmean([h1,h2]);
   
    xind=5; %location index for predictor 1(an example)
        x2ind=6; %location index for predictor 2
          y1=t0(xind);%location for 1st predictor
          y2=t00(x2ind);%location for predictor 2
          options = optimset('MaxFunEvals',30000,'MaxIter',90000);
            h=hh/sqrt(test11(xind)*test12(x2ind)); %adaptive bandwidth parameter      
              strt=zeros(1,m);
              d=fminsearch(@(c)formwteddensityregressionfromMh(c,X,[Y Y1],[y1 y2],N,M,Phi,t,s1,h),strt,options);
       
            d1=d(1:m);
          gamEst = FormGammaFromC(d1,Phi);
          gamDot = gradient(gamEst,1/(N));
          m1=locallinearregbivariate(XXX,Y,Y1,[y1 y2]);
          m1=(m1-A)/(B-A);
          fp=normpdf(t,m1,s1);fp=fp/(sum(fp)/N);
          fn = interp1(t, fp, (t(end)-t(1)).*gamEst + t(1)).*gamDot;
          fn=fn/(sum(fn)/N);
          fn=fn/(B-A);
          ft=exp(-y2)*(exp(-abs(t1-((2*y1 -1)^2))/(1/sqrt(2)))/(sqrt(2))) + (1-exp(-y2))*(normpdf(t1,(2*y1 +3),0.5));