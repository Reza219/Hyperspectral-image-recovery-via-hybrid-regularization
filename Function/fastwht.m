function Y = fastwht(X,perms,picks,N,M,mode)
% fast Walsh-Hadamard transform

nc=size(X,2);

switch mode
    case 1
        Y=zeros(M,nc);
        for ii=1:nc
            x = X(:,ii);
            fx = fWHtrans(x(perms))*sqrt(N);
            Y(:,ii) = fx(picks);
        end
    case 2
        Y=zeros(N,nc);
        for ii=1:nc
            fx = zeros(N,1);
            fx(picks) = X(:,ii);
            Y(perms,ii) = fWHtrans(fx)*sqrt(N);
        end
    otherwise
        error('Unknown mode');
end

end