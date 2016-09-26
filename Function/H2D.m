function Y=H2D(X,mode)
% 2-dimensional Haar transform

[nr,nc]=size(X);
nrr=sqrt(nr);
Y=zeros(nr,nc);

switch mode
    case 1
        for ii=1:nc
            T=haar_2d(reshape(X(:,ii),[nrr,nrr]));
            Y(:,ii)=T(:);
        end
    case 2
        for ii=1:nc
            T=haar_2d_inverse(reshape(X(:,ii),[nrr,nrr]));
            Y(:,ii)=T(:);
        end
end

end