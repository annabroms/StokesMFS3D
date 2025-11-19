function y = krylov_sqrt(A,z,m)


%Introduce a block method later? 
v(:,1) = z./norm(z);

for j =1:m
    w = A*v(:,j);
    if j>1
        w = w-H(j-1,j)*v(:,j-1);
    
    end
    H(j,j) = w'*v(:,j);
    if j<m
        w = w-H(j,j)*v(:,j);
        H(j+1,j) = norm(w);
        H(j,j+1) = norm(w);
        v(:,j+1) = w./H(j+1,j);
    end
end

e = zeros(m,1);
e(1) = 1; 
%make this function faster - avoid computing sqrt explicitly
y = norm(z)*v*sqrtm(H)*e;   


end