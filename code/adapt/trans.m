function Y2 = trans (Y, Ab, q)


N = size(Y,1); P = size(Y,2)/2;
blka = zeros(2*P,2*P);
b = Ab(4*P+1:end);
Y2 = zeros(N,2*P);
for p=1:P
    blka(2*(p-1)+1:2*p,2*(p-1)+1:2*p) = reshape(Ab(4*(p-1)+1: 4*p),2,2);
end
if q == 0
    G = zeros(2*P,2*P,2*P);
else
    G = randn(2*P,2*P,2*P)/(q);
end
for i=1:N
    z = Y(i,:)';
    for p=1:2*P
        Y2(i,p) = 0.5*z'*G(:,:,p)*z + blka(p,:)*z + b(p);
    end
end