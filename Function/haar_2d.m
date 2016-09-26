function Y=haar_2d(X)
% computes the Haar transform of an array

n=size(X,1);
s=1/sqrt(2);
Y=X;

k=n;
while 1<k
    k=floor(k/2);
    temp        =(Y(1:2:2*k-1,:)-Y(2:2:2*k,:))*s;
    Y(  1:  k,:)=(Y(1:2:2*k-1,:)+Y(2:2:2*k,:))*s;
    Y(k+1:2*k,:)=temp;
end

k=n;
while 1<k
    k=floor(k/2);
    temp        =(Y(:,1:2:2*k-1)-Y(:,2:2:2*k))*s;
    Y(:,  1:  k)=(Y(:,1:2:2*k-1)+Y(:,2:2:2*k))*s;
    Y(:,k+1:2*k)=temp;
end

end