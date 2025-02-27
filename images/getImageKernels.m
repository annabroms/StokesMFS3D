function Nimage = getImageKernels(singvec,rimage,nimage,rtest)
%% GETIMAGEKERNELS(rimage,nimage,rtest)

N1 = size(rimage,1); %number of source points
N2 = size(rtest,1); %Number of target points

if ~N1
    Nimage = [];
    return;
end

r = rimage'; 
rout = rtest'; 

%Pre-assign matrices
T = zeros(3*N1,3*N2); 
D = zeros(3*N1,3*N2); 

%parfor i = 1:N
if singvec(3) %Computation of the stresslet will be slow. This code can certainly be replaced. 
    %I don't do that now as we wont use the stresslet anyway.
    for i = 1:N1
        ind = 1:N2;
        z = r(:,i)-rout(:,ind);
        d = vecnorm(z,2,1);
        
        nvec = nimage(i,:);
        dotted = nvec*z;
    
        %Attempt at doing this faster, which seems to be much slower
        n = 3;
        j1 = N2;
        AA1 = [z; sparse(N2*3,j1)]; %pad with zeros and do sparse computation 
        AA1 = reshape(AA1,N2*3,[]);
        spAA1 = AA1(:,1:end-1);
        Atemp1 = spAA1*spAA1'; %this is rr'. A bit slow.. How can it be made faster? 
        
    
        linIndices1 = (0:n*((n*j1)+1):n*((n*j1)+1)*(j1-1))+reshape((1:n)'+n*j1*(0:n-1),[],1);
        %newA = reshape(A(linIndices),n,n,[]);
        A1 = reshape(Atemp1(linIndices1),n,j1*n);
        
        %% Build double layer (stresslet)
        C3 = 1./(d.^5);
        part2 = kron(C3.*dotted,ones(3)).*A1;
        
        %Mik = 1/8/pi*(part1+part2);
        
        DLP = 3/4/pi*(part2);
    
        T(3*(i-1)+1:3*i,1:end) = DLP;
    
        %% Build potential dipole
        C1 = 1./(d.^3);
        part1 = kron(C1,eye(3));
        C2 = 1./(d.^5);
        part2 = kron(C2,ones(3)).*A1;
        
        D(3*(i-1)+1:3*i,1:end) = (-part1+3*part2);
    
        
       
    end
end

%Build potential dipole

d1 = bsxfun(@minus,rimage(:,1)',rtest(:,1)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,rimage(:,2)',rtest(:,2));
d3 = bsxfun(@minus,rimage(:,3)',rtest(:,3));
r2 = d1.^2+d2.^2+d3.^2;   % dist^2 mat

ir = 1 ./ sqrt(r2);   % 1/r, incl wei
r4 = r2.*r2;
ir3 = ir./r2;
ir5 = ir ./ r4;        	% 1/r^3, incl wei


d1d2 = d1.*d2.*ir5;
d1d3 = d1.*d3.*ir5;
d2d3 = d2.*d3.*ir5;


D1 = zeros(3); D1(1,1) = 1;
D2 = zeros(3); D2(2,2) = 1;
D3 = zeros(3); D3(3,3) = 1;

E1 = zeros(3); E1(1,2) = 1; E1(2,1) = 1;
E2 = zeros(3); E2(1,3) = 1; E2(3,1) = 1;
E3 = zeros(3); E3(3,2) = 1; E3(2,3) = 1;

dxd_ir5 = kron(d1.^2.*ir5,D1) + kron(d2.^2.*ir5,D2) + kron(d3.^2.*ir5,D3) + kron(d1d2,E1) +...
    kron(d1d3,E2) + kron(d2d3,E3);
D = 1/4/pi*(kron(ir3,-eye(3)) + 3*dxd_ir5); %with scaling

%D = (kron(ir3,-eye(3)) + 3*dxd_ir5); 



%Could probably build the rotlet in a similar fashion... 

R2 = zeros(3*N2,3*N1);

for k = 1:N1
    r1= rotlet(rout',r(:,k)',[1 0 0]);
    r2= rotlet(rout',r(:,k)',[0 1 0]);
    r3= rotlet(rout',r(:,k)',[0 0 1]);
    r1 = r1';
    r2 = r2';
    r3 = r3';
    Rblock = [r1(:) r2(:) r3(:)];
    R2(:,(k-1)*3+1:3*k) = Rblock;
    % for l = 1:N2
    %     %Rblock = getRotletBlock(r(:,k)',rout(:,l)');
    %     Rblock = getRotletBlock(rout(:,l)',r(:,k)');
    %     R((k-1)*3+1:3*k,(l-1)*3+1:3*l) = Rblock;
    % 
    % end
end
R = 1/8/pi*R2;

%R = R2'; %not same correct scaling

%What singularities are to be included= sing_vec = S, R, T, D
Nimage = [];
if singvec(2)
    Nimage = R;
end
if singvec(3)
    Nimage = [Nimage T'];
end
if singvec(4)
    Nimage = [Nimage D];
end

% Nimage = [S' D'];
% %Nimage = [D'];
% %Nimage = [S' R'];
% Nimage = [R' S' D']; 
% Nimage = S';

%Nimage = [R' S']; 

end