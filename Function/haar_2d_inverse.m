function Y=haar_2d_inverse(X)
% computes the inverse Haar transform of an array

[n,m]=size(X);
s=1/sqrt(2);
Y=X;
T=zeros(n,m);

k=1;
while k*2<=n
    T(:,1:2:2*k-1)=(Y(:,1:k)+Y(:,1+k:2*k))*s;
    T(:,2:2:2*k)  =(Y(:,1:k)-Y(:,1+k:2*k))*s;
    Y(:,1:2*k)    =T(:,1:2*k);
    k=k*2;
end

k=1;
while k*2<=n
    T(1:2:2*k-1,:)=(Y(1:k,:)+Y(1+k:2*k,:))*s;
    T(2:2:2*k,:)  =(Y(1:k,:)-Y(1+k:2*k,:))*s;
    Y(1:2*k,:)    =T(1:2*k,:);
    k=k*2;
end

end