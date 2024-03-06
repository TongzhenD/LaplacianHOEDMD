function [U,Lambda,V] = hoedmd_stls(Psi,p,S)
%HOEDMD_STLS: the pth-order extended dynamic mode decomposition
%             based on the structured total least squares.
%Input:
%   Psi: stacked snapshots matrix, collects psi(F^{k-1}(L^n_t)), k=1,...,p+1
%   p: order of hoedmd
%   S: spectral complexity
%Output:
%   U: mode
%   Lambda: frequency
%   V: left eigenvectors
%by Weiyang Ding @Fudan January 25, 2024
M = round(size(Psi,1)/(p+1));% spatial dimension
indx = 1:(M*p);% x-part
indy = (M+1):(M*(p+1));% y-part
[P,~,~] = svd(Psi,'econ');
P = P(:,1:S);
[UP,Omega,VP] = svd(P(indx,:),'econ');
[U,Lambda,V] = eig(VP*(Omega\(UP'*P(indy,:))));
inprodUV = sum(conj(U).*V).*diag(Lambda)';
U = P((M+1):(M*2),:)*U;
V = UP*(Omega\(VP'*V));
normU = vecnorm(U);
U = bsxfun(@rdivide,U,normU);
V = bsxfun(@times,V,normU./inprodUV);
end

