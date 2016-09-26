function er = proposed_f(Im,L,np,Np,Ns,rp,rs,s,lam,gam1,gam2,ni,mss,N,perms1,perms2,picks11,picks22)

% problem size
Mp=round(rp*Np);
Ms=round(rs*Ns);
Mp1=round(.1*Np);
Ms1=round(.05*Ns);

% measurement matrix in spatial domain
picks1=picks11(1:Mp-Mp1);
Alp=@(Z,mode)lpfzigzag(Z,Np,np,Mp1,s,mode);
Arp=@(Z,mode)fastwht(Z,perms1,picks1,Np,Mp-Mp1,mode);
Ap=@(Z,mode)hybrid(Z,Alp,Arp,Mp1,Mp,mode);

% measurement matrix in spectral domain
picks2=picks22(1:Ms-Ms1);
Als=@(Z,mode)lpffwht(Z,Ns,Ms1,mode);
Ars=@(Z,mode)fastwht(Z,perms2,picks2,Ns,Ms-Ms1,mode);
As=@(Z,mode)hybrid(Z,Als,Ars,Ms1,Ms,mode);

% measurement
Y=Ap(As(Im,1)',1)'+N(1:Ms,1:Mp);

% initialization
X=Ap(As(Y,2)',2)';
Xn=X;
t=1;

% reconstruction
for ii=1:ni
    Xo=Xn;
    
    %l1 proximal
    X=L*wthresh(L'*X,'s',lam*gam2);
    
    %subgradient descent
    C=Y-Ap(As(X,1)',1)';
    X1=reshape(X',[np,np,Ns]);
    for jj=1:Ns
        X1(:,:,jj)=X1(:,:,jj)-lam*gam1*tv_l2_subgrad(X1(:,:,jj));
    end
    Xn=reshape(X1,[Np,Ns])'+lam*Ap(As(C,2)',2)';
    
    tk=t;
    t=(1+sqrt(1+4*tk^2))/2;
    X=Xn+(tk-1)/t*(Xn-Xo);
end

mse=mean(abs(X(:)-Im(:)).^2);
er=sqrt(mse/mss)*100;

end