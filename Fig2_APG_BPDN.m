%Fig 2, APG-BPDN

clear
addpath('Function')
load('.\DC.mat')
load('.\randsDC.mat')

Im=Im';
[L,~,~]=svd(Im(:,1:100:end),'econ');
[~,s,~]=jpgzzind(np,np);

%paramaters
lam=.25;                 %step-size
gam=.05;                 %regularization parameter
ni=200;                  %no. of iterations
rp=.1:.1:1;              %spatial measurement ratio
rs=[.1,.2,.5,1];         %spectral measurement ratio
mss=mean((Im(:)).^2);

figure
for rsi=rs
    er=[];
    rsi
    for rpi=rp
        rpi
        a=APG_BPDN_f(Im,L,np,Np,Ns,rpi,rsi,s,lam,gam,ni,mss,N,perms1,perms2,picks11,picks22);
        er=[er;a];
    end
    plot(rp*100,er)
    hold on
end