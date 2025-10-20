function tau_coarse = getCoarseSource(f,Rinv,Zi,Yi,R,opt)

db = opt.db;
N = opt.N;
M = opt.M;
P = opt.P;

proj_rhs = zeros(db*P,1);
for k = 1:P
    if opt.ellipsoid && (opt.lr>1)
        Rk = R{k};
        if opt.lr == 2
            res_vec = kron(eye(2),Rk)*Yi'*rotate_vector(f(3*(k-1)*M+1:3*k*M),Rk');
        else
            Yik = reshape(rotate_vector(Yi,Rk),3*M,db);
            res_vec = Yik'*f(3*(k-1)*M+1:3*k*M); 
        end
    else
        res_vec = Yi'*f(3*(k-1)*M+1:3*k*M);
    end
    proj_rhs(db*(k-1)+1:db*k) = res_vec;
end
coarse_coeff = Rinv*proj_rhs; % This matvec can be large if many basis functions
tau_coarse = zeros(N*P,1); 

for k = 1:P
   if opt.ellipsoid && (opt.lr>1)
       Rk = R{k};
       if opt.lr == 2
            res_vec = Zi*rotate_vector(coarse_coeff((k-1)*db+1:k*db),Rk');
            tau_coarse((k-1)*3*N+1:3*k*N) = rotate_vector(res_vec,Rk);
       else
            res_vec = reshape(rotate_vector(Zi,Rk),3*N,db)*coarse_coeff((k-1)*db+1:k*db);
            tau_coarse((k-1)*3*N+1:3*k*N) = res_vec; 
       end
       
   else
       tau_coarse((k-1)*3*N+1:3*k*N) = Zi*coarse_coeff((k-1)*db+1:k*db);
   end

end

end