function K = getKmat(r,q)
%GETKMAT(r,q) returns the matrix K that maps a vector of rigid body
%velocities of dimension 6P to a vector of velocities in the points r on
%the surfaces of the particles. For multiple particles, K is a block
%diagonal matrix where each block has 6 columns. 

cross_mat = @(x) -[0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
%for reference regarding the sign, see https://en.wikipedia.org/wiki/Cross_product


N = size(r,1); %
P = size(q,1); %
na = N/P; %number of points per particle

%starting at first molecule
m = 1;
debug = 1;


r = r'; 
q = q'; 

K = zeros(3*N,6*P);

%crete positive identity blocks
J = repmat(eye(3),na,1);


for i = 1:P %loop over particles
    B = zeros(3*na,3);
    %create particle block
    for k = 1:na
        B(3*(k-1)+1:3*k,:) = cross_mat(r(:,(i-1)*na+k)-q(:,i));    
    end
   
    %BB((i-1)*3*na+1:i*3*na,3*(i-1)+1:3*i) = B;
    K_block = [J,B];
    %from the python code: [3*offset:3*(offset+b.Nblobs), 6*k:6*k+6] 
    K((i-1)*3*na+1:i*3*na,6*(i-1)+1:6*i) = K_block;
end

end