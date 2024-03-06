function [U,Lambda,V] = hoedmd_gram(Psi0,Psi,p,S,ltra)
%HOEDMD_gram: An efficient implementation of Laplacian HOEDMD based on 
%             the Gram matrix 
%Input:
%   Psi0: the raw snapshots matrix, psi({L^n_t}), t=1,2,...T
%   Psi: stacked snapshots matrix, collects psi(F^{k-1}(L^n_t)), k=1,...,p+1
%   p: order of hoedmd
%   S: spectral complexity
%   ltra: list of lengths for each trajetory
%Output:
%   U: mode
%   Lambda: frequency
%   V: left eigenvectors
%by Tongzhen Dang @Fudan January 25, 2024
M=size(Psi0,1);
T=sum(ltra);
ind=zeros(1,size(Psi,2));
for k=1:length(ltra)
    c=sum(ltra(1:k));
    ind(c-ltra(k)-(k-1)*p+1:c-k*p)= 1+c-ltra(k):c-p;
end
G=Psi0'*Psi0; % Generate Gram matrix G
Gxx=zeros(size(ind,2));
Gxy=zeros(size(ind,2));
for j=1:p+1
    temp1=G(j:T-p+j-1,j:T-p+j-1);
    if j<p+1
        Gxx=Gxx+temp1(ind,ind); % Generate Gram matrix Gxx
        temp2=G(j:T-p+j-1,j+1:T-p+j);
        Gxy=Gxy+temp2(ind,ind); % Generate Gram matrix Gxy
    else
        Gp=temp1(ind,ind);
    end
end
indx = 1:(M*p);
[Q,sigma] = eig(Gxx+Gp,'vector');
[~,ind] = sort(sigma,'descend');
Q = Q(:,ind(1:S));
L = chol(Q'*Gxx*Q,'lower');
[U,Lambda,V] = eig(L'\(L\(Q'*Gxy*Q)));
inprodUV = sum(conj(U).*V).*diag(Lambda)';
U = Psi((M+1):(M*2),:)*(Q*U);
V = Psi(indx,:)*(Q*(L'\(L\V)));
normU = vecnorm(U);
U = bsxfun(@rdivide,U,normU);
V = bsxfun(@times,V,normU./inprodUV);
end

