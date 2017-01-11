clear all; 
close all;
N=100;
t=(1:N)/N;
n=100;%number of predictors

Y=rand(n,1);Y1=rand(n,1);% relevant predictors
k=18;% 18 nuisance/irrelevant predictors
Y2=rand(n,k);
Y11=exp(-Y1);
u=rand(n,1);
%X=((2*Y -1).^2+laprnd(n,1,0,1));
 X=((2*Y -1).^2+laprnd(n,1,0,1)).*(u<Y11) + (u>=Y11).*normrnd((3+2*Y),0.5);
% an example(bivariate)
XX=X;

A=min(X)- (std(X)/sqrt(n));
B=max(X)+ (std(X)/sqrt(n));
X=(X-A)/(B-A);

Z=cat(2,ones(n,1),Y);


T = 2*pi*(1:N)/N;
tt0=(1:N)/N;

M=ones(n,1);

%%% Nonparametric part starts here
%%fourier basis%%
K = 3; % number of predictors/2
Phi(1:K,:) = sqrt(2)*sin((1:K)'*T); % 3 sine components
Phi(K+1:2*K,:) = sqrt(2)*cos((1:K)'*T);% 3 cosine components
m = size(Phi,1);
options = optimset('MaxFunEvals',30000,'MaxIter',90000);

predictor=[Y Y1 Y2]; 


for i=1:n
   diff=predictor-ones(n,1)*predictor(i,:);
   normdiff=sqrt(sum(diff.^2,2));
   wt=normpdf(normdiff,0,1);
   wt=wt/sum(wt); %weight for nadaraya watson estimator
   M(i)=sum(wt.*XX); % the Nadaraya Watson estimate
   err(i)=XX(i)-M(i); % residual
   M(i)=(M(i)-A)/(B-A); % scaled estimate
end

s1=std(err)/(B-A); %scaled standard deviation estimate
figure(10); clf; 
for iter=1:10
y1=0.2;y2=rand(1,k);y=0.1;%location of predictors. Irrelevan predictors generated randomly
location=[y y1 y2];
for dim=1:(k+2) %k+2 because k irrelevant & 2 relevant predictors
    [test1(dim),test2(dim),hh(dim)]=ksdensity(predictor(:,dim),location(dim));
end
% the above finds the naive bandwidth estimate for each predictor. 
% It is computationally very fast but is not optimized through cross
% validation. Performance price paid for numerical efficiency
h1=harmmean(hh); %harmonic mean
h=h1/sqrt(prod(test1));%adaptive bandwidth at location
   strt=zeros(1,m);
    d=fminsearch(@(c)formwteddensityregressionfromMh(c,X,predictor,location,N,M,Phi,t,s1,h),strt,options);
   d1=d(1:m);
diff=predictor-ones(n,1)*location;
   normdiff=sqrt(sum(diff.^2,2));
   wt=normpdf(normdiff,0,1);
   wt=wt/sum(wt);
   m1=sum(wt.*XX);m1=(m1-A)/(B-A);%Nadaraya watson estimate at location
fp=normpdf(t,m1,s1);fp=fp/(sum(fp)/N);
gamEst = FormGammaFromC(d1,Phi);
gamDot = gradient(gamEst,1/(N));
fn = interp1(t, fp, (t(end)-t(1)).*gamEst + t(1)).*gamDot;
fn=fn/(sum(fn)/N);
fn=fn/(B-A);
fp=fp/(B-A);
t1=t.*(B-A) + A;%the grid for unscaled data

ft=exp(-y1)*(exp(-abs(t1-((2*y -1)^2))/(1/sqrt(2)))/(sqrt(2))) + (1-exp(-y1))*(normpdf(t1,(2*y +3),0.5));
%ft=exp(-abs(t1-((2*y -1)^2))/(1/sqrt(2)))/(sqrt(2));

plot(t1,ft,'k:',t1,fn,'k','LineWidth',2);
set(gca,'fontsize',18);
legend('True','Estimates','Location','northeast');
hold on;
end
hold off;
grid;
xlabel('Relevant predictor (0.1,0.2)')


