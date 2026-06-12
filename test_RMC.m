addpath PROPACK; % could disable the mex function inside

n = 2000; r = 5; 
alpha = 0.1; p = 40*r/n;

k = 2; % condition number
delta = (1-1/k)/(r-1);
Sig = diag(1:-delta:1/k);

ERR = zeros(2,10); % relative error in inf norm
ITER = zeros(2,10); % number of iterations

for rep = 1:10

    disp(rep)
    rng(rep)

    % define L
    P = randn(n,r);
    Q = randn(n,r);

    [U,~] = qr(P,'econ');
    [V,~] = qr(Q,'econ');

    % trim
    U = U.*repmat(min(ones(n,1),sqrt(r/n)./sqrt(sum(U.^2,2))),1,r);
    V = V.*repmat(min(ones(n,1),sqrt(r/n)./sqrt(sum(V.^2,2))),1,r);

    [U,~] = qr(U,'econ');
    [V,~] = qr(V,'econ');

    L = U*Sig*V';

    L_max = max(abs(L(:)));
    mu = max(max(sum(U.*U,2)),max(sum(V.*V,2)))*n/r;

    % define S
    IND = find(rand(n,n) < p);
    m = length(IND);

    idx = find(rand(m,1) < alpha);
    q = length(idx);
 
    mag = max(abs(L(:))); 
    s = 2*mag*rand(q,1)-mag;

    M = zeros(n,n);
    M(IND) = L(IND);
    M(IND(idx)) = M(IND(idx))+s;
    M = sparse(M); % converted to sparse matrix

    % RMC algorithms
    thresh = 1.1*mu*r/n; % should be bigger than L_max
    gamma = 0.7; % larger for harder problem

    [U1,Sig1,V1,~,~,err1] = RMC_SOFT(M,p,r,thresh,gamma);
    err1_inf = error_infty(U1*Sig1,V1',U*Sig,V')/L_max;

    ERR(1,rep) = err1_inf;
    ITER(1,rep) = length(err1);

    [U2,Sig2,V2,~,~,err2] = RMC_SCAD(M,p,r,thresh,gamma);
    err2_inf = error_infty(U2*Sig2,V2',U*Sig,V')/L_max;

    ERR(2,rep) = err2_inf;
    ITER(2,rep) = length(err2);

end

clearvars -except ERR ITER
