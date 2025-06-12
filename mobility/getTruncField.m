function u = getTruncField(P,q,rin_base,rout_base,lambda,Lcut)
%getTruncatedField(P,q,rin_base,rout_base,lambda,Lcut)

M = size(rout_base,1);
N = size(rin_base,1); 
u = zeros(P*M*3,1); 
vars.fmm = 0; 
for k = 1:P
    d = vecnorm(q(k,:)-q,2,2);
    d = d-2;
    ind = find(d<Lcut);

    %rin = rin_base+q(k,:);
    %lambda_in = lambda((k-1)*N*3+1:k*N*3);
    rout = rout_base+q(k,:);
   
    rin = [];
    lambda_in = [];
    for j = 1:length(ind)
        i = ind(j);
        rin = [rin; rin_base+q(i,:)];
        lambda_in = [lambda_in; lambda((i-1)*N*3+1:i*N*3)];
    end

    u((k-1)*3*M+1:k*3*M) = getFlow(lambda_in,rin,rout,vars);
       

end

end