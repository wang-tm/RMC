function [U,Sig,V,S_vec,si,err] = RMC_SCAD(M,p,r,thresh,gamma)
% inputs
% M: observed samples as a sparse matrix
% p: sampling ratio
% r: rank
% thresh: beta, initial threshold
% gamma: decay rate
%
% outputs
% U, Sig, V: rank-r svd of L
% S_vec: outlier estimate on observed locations
% si: success indicator, 1 if success
% err: relative error on observed locations in Frobenius norm
%
% Please refer to "Leave-One-Out Analysis for Nonconvex Robust Matrix 
% Completion with General Thresholding Functions" for details

[m,n] = size(M);
M_fro = norm(M,'fro');

[I,J,M_vec] = find(M);
q = length(M_vec);
IND = sub2ind([m n],I,J); % Omega

a = 3; % parameter for scad function, greater than 2

err = zeros(500,1);
err(1) = 1;

% initialization
S_vec = zeros(q,1);
idx1 = find(abs(M_vec) < 2*thresh);
S_vec(idx1) = sign(M_vec(idx1)).*max(abs(M_vec(idx1))-thresh,0);
idx2 = find(abs(M_vec) > a*thresh);
S_vec(idx2) = M_vec(idx2);
idx3 = setdiff((1:q)',[idx1;idx2]);
S_vec(idx3) = ((a-1)*M_vec(idx3)-a*thresh*sign(M_vec(idx3)))/(a-2);

D = sparse(I,J,(M_vec-S_vec)/p,m,n,q);
[U,Sig,V] = lansvd(D,r,'L'); % rank-r svd of a sparse matrix, factors of L^1

si = 0; % success indicator
for t = 1:499

    L_obs = compute_X_Omega(U*Sig,V,IND);
    D = sparse(I,J,L_obs+S_vec-M_vec,m,n,q);
    err(t+1) = norm(D,'fro')/M_fro;

    if err(t+1) < 1e-4
        si = 1;
        break
    else
        % update S
        D_vec = M_vec-L_obs;
        thresh_t = thresh*(gamma^t);
        S_vec = zeros(q,1);
        idx1 = find(abs(D_vec) < 2*thresh_t);
        S_vec(idx1) = sign(D_vec(idx1)).*max(abs(D_vec(idx1))-thresh_t,0);
        idx2 = find(abs(D_vec) > a*thresh_t);
        S_vec(idx2) = D_vec(idx2);
        idx3 = setdiff((1:q)',[idx1;idx2]);
        S_vec(idx3) = ((a-1)*D_vec(idx3)-a*thresh_t*sign(D_vec(idx3)))/(a-2);

        % update L
        D = sparse(I,J,(L_obs+S_vec-M_vec)/p,m,n,q);
        Vt = V'; Yforward = @(y) fast_multiply_left(U*Sig,Vt,D,y); clear Vt
        Ut = U'; Dt = D'; Ytranspose = @(y) fast_multiply_right(V*Sig,Ut,Dt,y); clear Ut Dt
        [U,Sig,V] = lansvd(Yforward,Ytranspose,m,n,r,'L');

    end

end

err = err(1:t+1);
