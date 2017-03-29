%Fig 2, proposed

clear
addpath('Function')

dataset='DC'; %'IP'; 'HI'; 'HO'; 'SD'; 'SF';

switch dataset
    case 'DC'
        load('.\DC.mat')
        load('.\randsDC.mat')
        gam1=.002;            %regularization parameter
        gam2=.004;            %regularization parameter
        
    case 'IP'
        load('.\IndianPines.mat')
        load('.\randsIP.mat')
        gam1=.001;
        gam2=.002;
        
    case 'HI'
        load('.\Harvard_i3.mat')
        load('.\randsHarv.mat')
        gam1=.001;
        gam2=.004;
        
    case 'HO'
        load('.\Harvard_oc4.mat')
        load('.\randsHarv.mat')
        gam1=.001;
        gam2=.004;
        
    case 'SD'
        load('.\StanfordDish.mat')
        load('.\randsStan.mat')
        gam1=.0002;
        gam2=.002;
        
    case 'SF'
        load('.\SanFrancisco.mat')
        load('.\randsStan.mat')
        gam1=.0002;
        gam2=.002;
end

Im=Im';
[L,~,~]=svd(Im(:,1:100:end),'econ');
[~,s,~]=jpgzzind(np,np);

%paramaters
lam=.25;                      %step-size
ni=150;                       %no. of iterations
rp=.1:.1:1;                   %spatial measurement ratio
rs=[.1,.2,.5,1];              %spectral measurement ratio
mss=mean((Im(:)).^2);

figure
for rsi=rs
    er=[];
    rsi
    for rpi=rp
        rpi
        a=proposed_f(Im,L,np,Np,Ns,rpi,rsi,s,lam,gam1,gam2,ni,mss,N,perms1,perms2,picks11,picks22);
        er=[er;a];
    end
    plot(rp*100,er)
    hold on
end