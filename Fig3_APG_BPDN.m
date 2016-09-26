%Fig 3, APG-BPDN

clear
addpath('Function')
load('.\DC.mat')
load('.\randsDC.mat')

Im=Im';
[L,~,~]=svd(Im(:,1:100:end),'econ');

% problem size
rp=.2;%.5;                    %spatial measurement ratio
rs=.1;%.2;                    %spectral measurement ratio
Mp=round(rp*Np);
Ms=round(rs*Ns);
Mp1=round(.1*Np);
Ms1=round(.05*Ns);

% measurement matrix in spatial domain
[~,s,~]=jpgzzind(np,np);
%perms1=randperm(Np);
%picks11=1+randperm(Np-1);
picks1=picks11(1:Mp-Mp1);
Alp=@(Z,mode)lpfzigzag(Z,Np,np,Mp1,s,mode);
Arp=@(Z,mode)fastwht(Z,perms1,picks1,Np,Mp-Mp1,mode);
Ap= @(Z,mode)hybrid(Z,Alp,Arp,Mp1,Mp,mode);

% measurement matrix in spectral domain
%perms2=randperm(Ns);
%picks22=1+randperm(Ns-1);
picks2=picks22(1:Ms-Ms1);
Als=@(Z,mode)lpffwht(Z,Ns,Ms1,mode);
Ars=@(Z,mode)fastwht(Z,perms2,picks2,Ns,Ms-Ms1,mode);
As= @(Z,mode)hybrid(Z,Als,Ars,Ms1,Ms,mode);

%noise
%N=.01*randn(Ns,Np);

% measurement
Y=Ap(As(Im,1)',1)'+N(1:Ms,1:Mp);

% paramaters
ni=194;%172;                  %no. of iterations
lam=.25;                      %step-size
gam=.05;                      %regularization parameter

% initialization
X=Ap(As(Y,2)',2)';
Xn=X;
t=1;
mse=zeros(ni,1);

% reconstruction
for ii=1:ni
    ii
    Xo=Xn;
    
    %l1 proximal
    X=L*wthresh(H2D((L'*X)',1)','s',lam*gam);
    X=H2D(X',2)';

    %gradient descent
    C=Y-Ap(As(X,1)',1)';
    Xn=X+lam*Ap(As(C,2)',2)';

    tk=t;
    t=(1+sqrt(1+4*tk^2))/2;
    X=Xn+(tk-1)/t*(Xn-Xo);
    
    mse(ii)=mean(abs(X(:)-Im(:)).^2);
end

mss=mean((Im(:)).^2);
er=sqrt(mse/mss)*100;
figure
plot(er)

X(X>1)=1;
X(X<0)=0;