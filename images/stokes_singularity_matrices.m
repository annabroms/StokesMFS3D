function [R,T,D] = stokes_singularity_matrices( ...
    sources,targets,stressletNormalDirections,activeOptionalFamilies,sourceNormals)
%STOKES_SINGULARITY_MATRICES Assemble dense free-space Stokes singularities.
%
% SOURCES is N-by-3 and TARGETS is M-by-3.  With unit viscosity and
% r = x-y, the implemented velocity kernels are
%
%   R(r;t)   = (t x r)/(8*pi*|r|^3),
%   T(r;n,q) = -6*r*(r.n)*(r.q)/(8*pi*|r|^5),
%   D(r;s)   = (-s/|r|^3 + 3*r*(r.s)/|r|^5)/(4*pi).
%
% ACTIVEOPTIONALFAMILIES is the logical row [useRotlets,useStresslets,
% usePotentialDipoles].  R and D use three interleaved coefficients per
% source.  T contains one 3N-column coefficient block per row of
% STRESSLETNORMALDIRECTIONS. If optional SOURCENORMALS is supplied as
% N-by-3, T instead contains one 3N-column block with an individual
% stresslet direction at each source. This avoids assembling three full
% Cartesian stresslet blocks merely to contract them afterward. Thus the
% output row order is
% [ux1;uy1;uz1;ux2;uy2;uz2;...], and R/D column order is
% [cx1;cy1;cz1;cx2;cy2;cz2;...].
%

validateattributes(sources,{'numeric'}, ...
    {'real','finite','2d','ncols',3},mfilename,'sources');
validateattributes(targets,{'numeric'}, ...
    {'real','finite','2d','ncols',3},mfilename,'targets');
validateattributes(stressletNormalDirections,{'numeric'}, ...
    {'real','finite','2d','ncols',3,'nonempty'},mfilename, ...
    'stressletNormalDirections');
if ~islogical(activeOptionalFamilies) || ...
        ~isequal(size(activeOptionalFamilies),[1,3])
    error('activeOptionalFamilies must be the logical row [useR,useT,useD].');
end

ns = size(sources,1);
nt = size(targets,1);
hasSrcNormals = nargin >= 5 && ~isempty(sourceNormals);
if hasSrcNormals
    validateattributes(sourceNormals,{'numeric'}, ...
        {'real','finite','2d','ncols',3,'nrows',ns}, ...
        mfilename,'sourceNormals');
    nd = 1;
else
    nd = size(stressletNormalDirections,1);
end
R = zeros(3*nt,3*ns*double(activeOptionalFamilies(1)));
T = zeros(3*nt,3*ns*nd*double(activeOptionalFamilies(2)));
D = zeros(3*nt,3*ns*double(activeOptionalFamilies(3)));
if ~any(activeOptionalFamilies)
    return
end

xr = 1:3:3*nt;
yr = 2:3:3*nt;
zr = 3:3:3*nt;
xc = 1:3:3*ns;
yc = 2:3:3*ns;
zc = 3:3:3*ns;

% The M-by-N pair grids are small compared with the requested dense
% 3M-by-3N blocks and avoid recomputing geometry for every vector component.
rx = targets(:,1)-sources(:,1).';
ry = targets(:,2)-sources(:,2).';
rz = targets(:,3)-sources(:,3).';
r2 = rx.^2+ry.^2+rz.^2;
if any(r2 == 0,'all')
    [it,is] = find(r2 == 0,1);
    error('stokes_singularity_matrices:CoincidentPoints', ...
        'Target %d coincides with source %d.',it,is);
end
ir3 = r2.^(-3/2); % |r|^{-3}

if activeOptionalFamilies(1)
    cR = ir3/(8*pi);
    R(xr,yc) = cR.*rz;
    R(xr,zc) = -cR.*ry;
    R(yr,xc) = -cR.*rz;
    R(yr,zc) = cR.*rx;
    R(zr,xc) = cR.*ry;
    R(zr,yc) = -cR.*rx;
end

if activeOptionalFamilies(2) || activeOptionalFamilies(3)
    ir5 = ir3./r2; % |r|^{-5}
end

if activeOptionalFamilies(2)
    for d = 1:nd
        if hasSrcNormals
            rDotN = rx.*sourceNormals(:,1).'+ ...
                ry.*sourceNormals(:,2).'+rz.*sourceNormals(:,3).';
        else
            n = stressletNormalDirections(d,:);
            rDotN = n(1)*rx+n(2)*ry+n(3)*rz;
        end
        cT = (-6/(8*pi))*rDotN.*ir5;
        c0 = (d-1)*3*ns;
        T(xr,c0+xc) = cT.*rx.*rx;
        T(xr,c0+yc) = cT.*rx.*ry;
        T(xr,c0+zc) = cT.*rx.*rz;
        T(yr,c0+xc) = cT.*ry.*rx;
        T(yr,c0+yc) = cT.*ry.*ry;
        T(yr,c0+zc) = cT.*ry.*rz;
        T(zr,c0+xc) = cT.*rz.*rx;
        T(zr,c0+yc) = cT.*rz.*ry;
        T(zr,c0+zc) = cT.*rz.*rz;
    end
end

if activeOptionalFamilies(3)
    diagD = -ir3/(4*pi);
    outerD = 3*ir5/(4*pi);
    D(xr,xc) = outerD.*rx.*rx+diagD;
    D(xr,yc) = outerD.*rx.*ry;
    D(xr,zc) = outerD.*rx.*rz;
    D(yr,xc) = outerD.*ry.*rx;
    D(yr,yc) = outerD.*ry.*ry+diagD;
    D(yr,zc) = outerD.*ry.*rz;
    D(zr,xc) = outerD.*rz.*rx;
    D(zr,yc) = outerD.*rz.*ry;
    D(zr,zc) = outerD.*rz.*rz+diagD;
end
end
