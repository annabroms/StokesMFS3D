function M = generate_RPY_matrix(r,a)
%%GENERATE_RPY_MATRIX(r,q,a) generates the blob-blob matrix N for all
%%"blobs" in the system, where the blob radius (RPY regularisation parameter) is a and the blobs positions
%%are given in the matrix r of size (number of blobs x 3). The matrix q
%%contains the center coordinate for all the multiblob particles. 

%Create main block of the mobility: 

N = size(r,1); 
r = r'; 
%MM = 1/6/pi/a*eye(3*N); %is this correct? 
MM = 1/12/pi/a*eye(3*N); %is this correct? 

%can this be done with e.g. a parfor loop instead and then reshaping to a
%matrix?
for i = 1:N
    
    %construct upper part of system matrix
    for k = i+1:N
        
        Mik = mobility_RPY_block(r(:,i)-r(:,k),a); 
        %particle 1-particle 2 is the convention for the RPY tensor
            
        %MM(3*(i-1)+1:3*i,3*(k-1)+1:3*k) = Mik;
        %MM(3*(k-1)+1:3*k,3*(i-1)+1:3*i) = Mik;  
        
        MM(3*(i-1)+1:3*i,3*(k-1)+1:3*k) = Mik;
        
    end     
    
end
M = (MM+MM');


end

function  M12 = mobility_RPY_block(r,a)
%r = (2+gap)*a*d;
%should use a as a parameter here
if nargin == 1
    a = 1;
end

%should use a as a parameter here
% C1 = @(x) 3*a./(4*x)+a^3./(2*x.^3);
% C2 = @(x) 3*a./(4*x)-3*a^3./(2*x.^3);
% C3 = @(x) 1-9*x/(32*a);
% C4 = @(x) 3*x/(32*a);

% r = reshape(r,3,1);
% A = r*r';
A = r(:)*r(:)';
d = norm(r);

if d > 2*a
    %M12 = 1/(6*pi*a)*(C1(d)*eye(3)+C2(d)*A./d^2);
    %M12 = 1/(6*pi*a)*(C1(d,a)*eye(3)+C2(d,a)*A./d^2);
    M12 = 1/(6*pi*a)*((3*a./(4*d)+a^3./(2*d.^3))*eye(3)+...
        (3*a./(4*d)-3*a^3./(2*d.^3))*A./d^2);
else
    %M12 = 1/(6*pi*a)*(C3(d)*eye(3)+C4(d)*A./d^2);
    M12 = 1/(6*pi*a)*(C3(d,a)*eye(3)+C4(d,a)*A./d^2);
end



end

function res = C1(x,a)
    res = 3*a./(4*x)+a^3./(2*x.^3);
end

function res = C2(x,a)
    res = 3*a./(4*x)-3*a^3./(2*x.^3);
end

function res = C3(x,a)
    res = 1-9*x/(32*a);
end

function res = C4(x,a)
    res = 3*x/(32*a);
end
