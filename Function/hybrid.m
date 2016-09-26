function Y = hybrid(X,Al,Ar,M1,M,mode)
% constructs hybrid measurement matrix

switch mode
    case 1
        Y=[Al(X,1);Ar(X,1)];
    case 2
        Y=Al(X(1:M1,:),2)+Ar(X(M1+1:M,:),2);
    otherwise
        error('Unknown mode');
end

end