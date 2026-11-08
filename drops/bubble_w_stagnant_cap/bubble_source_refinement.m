function patch = bubble_source_refinement(geometry,opt,validationLevel)
%BUBBLE_SOURCE_REFINEMENT Sample local source footprints on both edge sides.
% validationLevel=0 gives fit points; positive levels give independent audits.
% For depth d, meridional offsets are radialFactors*d and azimuthal arc offsets
% are azimuthFactors*d. Sources and the existing global grid are unchanged.
if nargin<3, validationLevel = 0; end
radial = opt.sourceRefinementRadialFactors;
azimuth = opt.sourceRefinementAzimuthFactors;
if validationLevel>0
    % Denser audits include new distances and the exact source azimuths.
    radial = logspace(log10(min(radial)/2),log10(max(radial)*2), ...
        2*validationLevel*numel(radial)+1);
    azimuth = unique([linspace(min(azimuth),max(azimuth), ...
        2*validationLevel*max(1,numel(azimuth)-1)+1),0]);
end
firstLevel = max(1,opt.n_level-opt.sourceRefinementLevels+1);
sources = find(geometry.sourceRingIndex>=firstLevel);
maxCount = numel(sources)*2*numel(radial)*numel(azimuth);
points = zeros(maxCount,3);
sourceIndices = zeros(maxCount,1);
depths = zeros(maxCount,1);
count = 0;
for index = sources(:).'
    y = geometry.sourcePoints(index,:)-opt.centre;
    depth = opt.radius-norm(y);
    psi = atan2(y(2),y(1));
    for factor = radial
        for side = [-1 1]
            theta = opt.phi+side*factor*depth/opt.radius;
            if theta<=0 || theta>=pi, continue; end
            angles = psi+azimuth*depth/(opt.radius*sin(theta));
            selected = count+(1:numel(angles));
            unit = [sin(theta)*cos(angles(:)),sin(theta)*sin(angles(:)), ...
                cos(theta)*ones(numel(angles),1)];
            points(selected,:) = opt.centre+opt.radius*unit;
            sourceIndices(selected) = index;
            depths(selected) = depth;
            count = count+numel(angles);
        end
    end
end
points = points(1:count,:);
% Preserve every old node and remove repeated candidates before appending.
tolerance = 1e-11*opt.radius;
[~,indices] = uniquetol(points,tolerance,'ByRows',true,'DataScale',1);
indices = sort(indices);
isOld = ismembertol(points(indices,:),geometry.points,tolerance, ...
    'ByRows',true,'DataScale',1);
indices = indices(~isOld);
patch.points = points(indices,:);
patch.normals = (patch.points-opt.centre)/opt.radius;
patch.theta = acos(max(-1,min(1,patch.normals(:,3))));
patch.cap = patch.theta<opt.phi;
patch.mobile = ~patch.cap;
patch.edgeDistance = opt.radius*abs(patch.theta-opt.phi);
patch.sourceIndex = sourceIndices(indices);
patch.sourceDepth = depths(indices);
patch.count = numel(indices);
assert(all(patch.edgeDistance>tolerance), ...
    'bubble_source_refinement:UnresolvedPatch', ...
    'Source-depth patches are too close to the edge for floating-point geometry.');
end
