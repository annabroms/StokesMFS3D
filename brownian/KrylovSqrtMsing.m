function [y,iter_errv, exact_errv] = KrylovSqrtMsing(M,z,tol,yexact)
%Borrowed from FBIM

if isa(M,'function_handle')
   afun = M; 
else
   afun = @(x) M*x;
end

iter_err = 1;
v = z/norm(z);
y = v;

%yexact = sqrtm((M+M')/2)*z;
%yexact = real(yexact);

j=1;
while iter_err > tol
    
    w = afun(v(:,j));
    
    h(j,j) = w'*v(:,j);
    
    
    if j>1
       w = w - h(j-1,j)*v(:,j-1); 
    end
    
    w = w - h(j,j)*v(:,j);
    
    h(j+1,j) = norm(w);
    h(j,j+1) = h(j+1,j);
    
    v(:,j+1) = w/norm(w);
    
    hsqrt = sqrtm(h(1:j,1:j));
    
    yold = y;
    y = norm(z) * ( v(:,1:j)*hsqrt(:,1));
    y = real(y);
    
    iter_err = norm(y-yold)/norm(yold);
    
    iter_errv(j) = iter_err;
    
    if nargin == 4
        exact_errv(j) = norm(y-yexact)/norm(yexact);
    end
    
    j = j+1;
%    [j,iter_err]    
end

if nargin<4
    exact_errv = [];
end

end