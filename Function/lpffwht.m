function Y = lpffwht(X,N,M,mode)
% low-pass filtering using fast Walsh-Hadamard transform

nc=size(X,2);

switch mode
    case 1
        Y=zeros(M,nc);
        for ii=1:nc
            fx = fWHtrans(X(:,ii));
            Y(:,ii) = fx(1:M)*sqrt(N);
        end
    case 2
        Y=zeros(N,nc);
        for ii=1:nc
            x = [X(1:M,ii);zeros(N-M,1)];
            Y(:,ii) = fWHtrans(x)*sqrt(N);
        end
    otherwise
        error('Unknown mode');
end

end