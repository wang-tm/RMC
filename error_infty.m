function err = error_infty(U1,V1t,U2,V2t)

n = size(V1t,2);

err = zeros(1,n);
for i = 1:n
    v1 = U1*V1t(:,i);
    v2 = U2*V2t(:,i);
    err(i) = max(abs(v1-v2));
end

err = max(err);