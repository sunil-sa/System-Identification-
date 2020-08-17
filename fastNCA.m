function [A, P] = fastNCA(Z, Astruct, p);
[nsamples, nvar] = size(Z);
[u s v] = svd(Z,'econ');
W = u(:,1:p);
Amix = zeros(nsamples,p);
for k = 1:p
    [Zc Zr] = rearrange(W,Astruct,k);
    [u s v] = svd(Zr);
    S = v(:,p);
    Zcp= Zc*S;
    [u s v] = svd(Zcp,'econ');
    ak = u(:,1);
    Amix(1:size(Zc,1),k) = ak;
end
A = reconstitute(Amix,Astruct);
P = inv(A'*A)*A'*Z;