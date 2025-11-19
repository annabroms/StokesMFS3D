function T = getTraction(r, rout, nn)
%GETTRACTION  Evaluate 3D traction matrix for Stokes single-layer kernel.
%
%   T = GETTRACTION(r, rout, nn)
%
%   Constructs the 3×3N traction matrix mapping source forces at locations `r`
%   to surface tractions (stress dotted with normal) at target points `rout`
%   with outward unit normals `nn`. This corresponds to applying the 
%   Stokeslet traction kernel in 3D free-space.
%
%   INPUTS:
%     r      - N×3 array of source points (e.g., inner surface nodes).
%     rout   - M×3 array of target points (e.g., outer surface nodes).
%     nn     - M×3 array of unit normals at each target point in `rout`.
%
%   OUTPUT:
%     T      - 3M×3N traction matrix, such that T * f approximates the 
%              surface traction field induced by source forces f.
%
%   NOTES:
%     - This routine assumes Stokes flow in 3D, with viscosity μ = 1.
%     - The matrix maps from 3N source force densities to 3M traction values.
%     - The implementation uses the singular fundamental solution and assumes
%       sources and targets are distinct (no regularization).
%
%   See also: getStokeslet, getFlow, getStresslet
%
% Anna Broms, Nov 19, 2025


%There is fmm code for doing this. Try that!

N1 = size(r,1); %number of source points
N2 = size(rout,1); %Number of target points

r = r'; 
rout = rout'; 

%Pre-assign matrices
T = zeros(3*N1,3*N2); 

%
for i = 1:N1 %Loop over sources 
    ind = 1:N2; % Take out all targets
    z = r(:,i)-rout(:,ind);
    d = vecnorm(z,2,1);
    
    %nvec = nn(i,:);
    nvec = nn(ind,:);
   % dotted = nvec*z;
    dotted = sum(nvec'.*z);

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
    %DLP = -6*(part2); %This can be compared to what Alex is doing for FMM I think
    %DLP = 
    %DLP = 3/8/pi*(part2); %only for debug

    T(3*(i-1)+1:3*i,1:end) = DLP;
end

%% Cleaner code that is producing the same thing: 
% Nsrc = numel(xsrc);
% Ntar = numel(xtar);
% T = zeros(3*Ntar, 3*Nsrc);
% 
% prefac = -6 / (8*pi*mu);
% 
% for j = 1:Nsrc
%     Vector differences r = x - y_j
%     rx = xtar - xsrc(j);
%     ry = ytar - ysrc(j);
%     rz = ztar - zsrc(j);
%     r2 = rx.^2 + ry.^2 + rz.^2;
%     r5 = r2.^(2.5);
% 
%     Dot with normal at target
%     rdotn = rx.*nx + ry.*ny + rz.*nz;
% 
%     Compute 3x3 block per target
%     Txx = prefac * rdotn .* (rx.*rx) ./ r5;
%     Txy = prefac * rdotn .* (rx.*ry) ./ r5;
%     Txz = prefac * rdotn .* (rx.*rz) ./ r5;
%     Tyx = prefac * rdotn .* (ry.*rx) ./ r5;
%     Tyy = prefac * rdotn .* (ry.*ry) ./ r5;
%     Tyz = prefac * rdotn .* (ry.*rz) ./ r5;
%     Tzx = prefac * rdotn .* (rz.*rx) ./ r5;
%     Tzy = prefac * rdotn .* (rz.*ry) ./ r5;
%     Tzz = prefac * rdotn .* (rz.*rz) ./ r5;
% 
%     %Assemble into block structure
%     for i = 1:Ntar
%         T(3*i-2:3*i, 3*j-2:3*j) = [ ...
%             Txx(i), Txy(i), Txz(i);
%             Tyx(i), Tyy(i), Tyz(i);
%             Tzx(i), Tzy(i), Tzz(i)];
%     end
% end

     
end
