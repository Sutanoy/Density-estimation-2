clear all; 
close all;


N=100;% number of grid points
t=(1:N)/N; % grid points for scaled responses(scaled to [0,1]
n=100;% number of observations

T = 2*pi*(1:N)/N;%evaluation points for the fourir basis elements


K = 3;%(Number of basis elements/2). We have used ad hoc 6 basis elements for all analysis.
% One can obtain a desired number of basis elements using Algorithm 1 as described for densiy
% estimation. But that stage was skipped for ease of illustration of the
% main ideas.


Phi(1:K,:) = sqrt(2)*sin((1:K)'*T); %sine components of Fourier basis elements
Phi(K+1:2*K,:) = sqrt(2)*cos((1:K)'*T); %cosine components of the same
m = size(Phi,1);
t0=0:0.1:1; % a grid over the predictor support on which we performed analysis for comparative purposes.

err=zeros(n,1);
   
   Y=rand(n,1); %though not required, all simulated examples used predictors generated from uniform density estimation.
   X=(2*Y -1).^2+laprnd(n,1,0,1); % an example of simulated data, from Laplace distribution
   XXX=X;
   A=min(X)- (std(X)/sqrt(n)); % estimation of bounds(lower) 
   B=max(X)+ (std(X)/sqrt(n)); %estimation of bounds(upper)
   X=(X-A)/(B-A);%scaling
   t1=t.*(B-A) + A;%grid for unscaled data
   
  M=zeros(n,1);
   for i=1:n
       M(i)=locallinearreg(XXX,Y,Y(i,1)); 
       % This portion is for obtaining the initial density mean function. 
       %Often this step gets slow and one can replace this step with
       %Nadaraya watson or any other naively chosen estimate
       err(i)=XXX(i)-M(i); %residual estimates
       M(i)=(M(i)-A)/(B-A); %scaling the men function for scaled data
   end
   
   s1=std(err)/(B-A); %scaling the standard deviation estimate
   [test1,test2,hh]=ksdensity(Y,t0);
   %The initial naive estimate of the smoothness of the data
   xind=5;
   %an illustrative and arbitrary location of the predictor 
   %where we perform the conditional density estimation
       y=t0(xind); %the location of the predictor Y
       options = optimset('MaxFunEvals',30000,'MaxIter',90000);

       h=hh/sqrt(test1(xind)); 
       % the adaptive bandwidth parameter at the location
   m1=locallinearreg(XXX,Y,y); % the initial density mean function
       m1=(m1-A)/(B-A);
       
        strt=zeros(1,m); %starting point, origin of tangent space
        d=fminsearch(@(c)formwteddensityregressionfromMh(c,X,Y,y,N,M,Phi,t,s1,h),strt,options);
       %the function formwteddensityregressionfromMh computes the objective
       %function for different values of c. Here c refers to the
       %coefficients of the basis elements in the tangent space
       %representation
   d1=d(1:m);% the optimal values as obtained through FMINSEARCH
   gamEst = FormGammaFromC(d1,Phi);%obtain the corresponding \gamma from the tangent space representation
   gamDot = gradient(gamEst,1/(N));%obtain \dot{\gamma} from \gamma
       
       fp=normpdf(t,m1,s1);fp=fp/(sum(fp)/N);%initial conditional density estimate shape
       fn = interp1(t, fp, (t(end)-t(1)).*gamEst + t(1)).*gamDot;%warping process 
       fn=fn/(sum(fn)/N);%normalizing for numerical accuracy
       fn=fn/(B-A);%scaling back to UNscaled data
       fp=fp/(B-A);%same for ease of comparison with the warped shape
   ft=exp(-abs(t1-((2*y -1)^2))/(1/sqrt(2)))/(sqrt(2)); % the corresponding illustrative true density shape
   figure(1); clf; 
plot(t1, fp, 'b', t1, fn,'r', t1, ft, 'g', 'LineWidth',3); 
set(gca,'fontsize',18);
grid; 
legend('initial','warped','true','Location','northeast'); 

