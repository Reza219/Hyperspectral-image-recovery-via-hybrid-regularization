function Y = lpfzigzag(X,Np,np,Mp,s,mode)
% low-pass filtering using zigzag pattern and fast Walsh-Hadamard transform

switch mode
    case 1
        n2=size(X,2);
        Y=zeros(Mp,n2);
        for ii=1:n2
            Xi=reshape(X(:,ii),[np,np]);
            T1=zeros(np,np);
            T2=zeros(np,np);
            for jj=1:np
                T1(:,jj)=fWHtrans(Xi(:,jj));
            end
            for jj=1:np
                T2(jj,:)=fWHtrans(T1(jj,:));
            end
            for jj=1:Mp
                Y(jj,ii)=T2(s(jj,1),s(jj,2))*np;
            end
        end
    case 2
        [n1,n2]=size(X);
        Y=zeros(Np,n2);
        for ii=1:n2
            T1=zeros(np,np);
            T2=zeros(np,np);
            T3=zeros(np,np);
            for jj=1:n1
                T1(s(jj,1),s(jj,2))=X(jj,ii);
            end
            for jj=1:np
                T2(:,jj)=fWHtrans(T1(:,jj));
            end
            for jj=1:np
                T3(jj,:)=fWHtrans(T2(jj,:));
            end
            Y(:,ii)=T3(:)*np;
        end
    otherwise
        error('Unknown mode');
end

end