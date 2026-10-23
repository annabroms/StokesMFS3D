%OCT17_TWO_SPHERE_IMAGE_FIELD_COMPARISON Pressure and traction at delta=0.1.
% Run this script deliberately: the refined one-body SVDs are substantial.
% All cases have the same geometry, rigid velocities, viscosity, and zero
% background flow. The source vectors stored in details.lambda are the
% unit-viscosity velocity coefficients; image_pressure/image_traction return
% physical fields at the viscosity below.

clear;
close all;
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root,'setpath.m'));

mu = 1.7;
q = [0,0,0;2.1,0,0];
U = [1;0.4;0;0;0;0.25; -1;-0.4;0;0;0;-0.25];
% imageFamilies tokens: S Stokeslet, R rotlet, D potential dipole, T legacy
% stresslet, T_t line-parallel, T_n first transverse, T_n2 second transverse
% stresslet, Q_s trace-free axial Stokeslet quadrupole. T excludes T_*.
imageFamilies = "SRD";
imageCounts = [3,6,10];       % Eq. (39) selects 3 at this gap.
svdCutoff = 5e-12;
tightTolerance = 1e-9;
looseTolerance = 1e-6;

base = struct('viscosity',mu,'collocationRatio',1.2, ...
    'svdCutoff',svdCutoff,'gmresTolerance',tightTolerance, ...
    'maxIterations',1000,'solverMode',"oneBody", ...
    'validationCount',2000,'validationGapCount',1001, ...
    'batchSize',128);
referenceCounts = [1353,1923,2739,3529];
cases = repmat(struct('name',"",'opt',struct(),'F',[], ...
    'it',0,'lambdaNorm',0,'boundaryResidual',0,'details',struct(), ...
    'elapsed',0),numel(referenceCounts)+numel(imageCounts)+1,1);
nRef = numel(referenceCounts);
refIndex = nRef;
baseImageIndex = nRef+1;
lastImageIndex = nRef+numel(imageCounts);
looseIndex = numel(cases);
for j = 1:nRef
    opt = base;
    opt.useImages = false;
    opt.N_image = 0;
    opt.proxyRadius = 0.75;
    opt.N = referenceCounts(j);
    cases(j).name = "reference N="+string(opt.N);
    cases(j).opt = opt;
end
for j = 1:numel(imageCounts)
    opt = base;
    opt.N = 686;
    opt.proxyRadius = 0.63;
    opt.useImages = true;
    opt.imageFamilies = imageFamilies; % See token legend above.
    opt.N_image = imageCounts(j);
    opt.capMode = "fibonacci";
    opt.capAngles = [pi/5,pi/60];
    opt.capUpsampling = [6,6];
    opt.rowWeightMode = "paper";
    index = nRef+j;
    cases(index).name = "images "+imageFamilies+ ...
        ", Nimage="+string(imageCounts(j));
    cases(index).opt = opt;
end
cases(looseIndex).name = "images "+imageFamilies+", Nimage="+ ...
    string(imageCounts(1))+", loose GMRES";
cases(looseIndex).opt = cases(baseImageIndex).opt;
cases(looseIndex).opt.gmresTolerance = looseTolerance;

fprintf('Two-sphere field comparison: gap %.3g, viscosity %.3g\n', ...
    norm(q(2,:)-q(1,:))-2,mu);
fprintf('Image accumulation radius %.12f, reference proxy radius %.3f\n', ...
    1+0.1/2-sqrt(0.1+0.1^2/4),cases(1).opt.proxyRadius);
for j = 1:numel(cases)
    tic;
    [cases(j).F,cases(j).it,cases(j).lambdaNorm, ...
        cases(j).boundaryResidual,cases(j).details] = ...
        solve_resistance_with_images(q,U,cases(j).opt);
    cases(j).elapsed = toc;
    body = cases(j).details.grid.body;
    sourceCounts = arrayfun(@(b) sum(arrayfun( ...
        @(z) size(z.points,1),b.basis)),body);
    collocationCounts = arrayfun(@(b) size(b.collocation,1),body);
    fprintf(['%s: sources [%s], collocation [%s], system %d x %d, ', ...
        'SVD cutoff %.1e, ranks [%s], GMRES requested %.1e, ', ...
        'achieved %.3e, it %d, boundary %.3e, |lambda| %.3e, ', ...
        'time %.1f s\n'],cases(j).name, ...
        num2str(sourceCounts(:).'),num2str(collocationCounts(:).'), ...
        cases(j).details.systemSize,svdCutoff, ...
        num2str(cases(j).details.retainedRanks.'), ...
        cases(j).opt.gmresTolerance, ...
        cases(j).details.trueGmresResidual,cases(j).it, ...
        cases(j).boundaryResidual,cases(j).lambdaNorm,cases(j).elapsed);
end

offsets = [0.005,0.025,0.05,0.1,0.5,2];
gridRequests = [2000,5000];
capCounts = [1001,2001];
field = cell(2,1);
targetSet = cell(2,1);
for level = 1:2
    [design,~] = get_sphdesign(gridRequests(level));
    surface = zeros(0,3);
    normals = zeros(0,3);
    desired = zeros(0,3);
    for body = 1:2
        direction = (q(3-body,:)-q(body,:))/norm(q(3-body,:)-q(body,:));
        sphere = struct('axes',ones(3,1),'Q',eye(3),'c',q(body,:).');
        capOpt = struct('collocationPatchMode',"fibonacci", ...
            'collocationPatchPointCount',capCounts(level), ...
            'collocationPatchAngle',pi/5);
        [cap,~] = collocation_surface_patch( ...
            sphere,q(body,:).'+direction.',capOpt);
        points = [q(body,:)+design;cap];
        nn = points-q(body,:);
        surface = [surface;points]; %#ok<AGROW>
        normals = [normals;nn]; %#ok<AGROW>
        desired = [desired;reshape( ...
            getKmat(points,q(body,:))*U(6*(body-1)+(1:6)),3,[]).']; %#ok<AGROW>
    end
    targetSet{level}.surface = surface;
    targetSet{level}.normals = normals;
    targetSet{level}.desired = desired;
    targetSet{level}.designCount = size(design,1);
    targetSet{level}.capCount = capCounts(level);

    exterior = zeros(0,3);
    labels = strings(0,1);
    for distance = offsets
        points = surface+distance*normals;
        valid = vecnorm(points-q(1,:),2,2)>1+1e-12 & ...
            vecnorm(points-q(2,:),2,2)>1+1e-12;
        exterior = [exterior;points(valid,:)]; %#ok<AGROW>
        labels = [labels;repmat("offset "+string(distance), ...
            nnz(valid),1)]; %#ok<AGROW>
    end
    axisPoints = [1.005,1.025,1.05,1.075,1.095];
    transverse = linspace(-0.35,0.35,21);
    [yy,zz] = ndgrid(transverse,transverse);
    for axisPoint = axisPoints
        points = [axisPoint*ones(numel(yy),1),yy(:),zz(:)];
        valid = vecnorm(points-q(1,:),2,2)>1+1e-12 & ...
            vecnorm(points-q(2,:),2,2)>1+1e-12;
        exterior = [exterior;points(valid,:)]; 
        labels = [labels;repmat("gap",nnz(valid),1)]; 
    end
    targetSet{level}.exterior = exterior;
    targetSet{level}.labels = labels;

    if level == 1
        selected = 1:numel(cases);
    else
        selected = [refIndex,baseImageIndex,lastImageIndex];
    end
    field{level} = repmat(struct('p',[],'t',[],'uResidual',NaN), ...
        numel(cases),1);
    for j = selected
        details = cases(j).details;
        sourceGrid = details.grid;
        lambda = details.lambda;
        field{level}(j).p = image_pressure(sourceGrid,lambda,exterior,mu);
        field{level}(j).t = image_traction(sourceGrid,lambda, ...
            surface,normals,mu);
        velocity = zeros(size(surface));
        for body = 1:2
            basis = sourceGrid.body(body).basis;
            coeff = image_source_coefficients(basis,lambda{body});
            for first = 1:128:size(surface,1)
                rows = first:min(first+127,size(surface,1));
                velocity(rows,:) = velocity(rows,:)+ ...
                    direct_image_velocity(basis,coeff,surface(rows,:));
            end
        end
        field{level}(j).uResidual = max(vecnorm( ...
            velocity-desired,2,2))/max(vecnorm(desired,2,2));
    end
    fprintf(['Validation level %d: %d design + %d cap points/body, ', ...
        '%d exterior pressure points, %d traction points\n'], ...
        level,size(design,1),capCounts(level),size(exterior,1), ...
        size(surface,1));
end

for level = 1:2
    ref = field{level}(refIndex);
    pScale = max(abs(ref.p));
    tScale = max(vecnorm(ref.t,2,2));
    fprintf(['Validation level %d: pressure scale %.3e, traction scale ', ...
        '%.3e, reference independent BC %.3e\n'], ...
        level,pScale,tScale,ref.uResidual);
    if level == 1
        selected = [1:refIndex-1,baseImageIndex:numel(cases)];
    else
        selected = [baseImageIndex,lastImageIndex];
    end
    for j = selected
        dp = abs(field{level}(j).p-ref.p);
        dt = vecnorm(field{level}(j).t-ref.t,2,2);
        fprintf(['%s: max |dp| %.3e, global relative %.3e; ', ...
            'max |dt| %.3e, global relative %.3e; ', ...
            'independent BC %.3e\n'],cases(j).name, ...
            max(dp),max(dp)/pScale,max(dt),max(dt)/tScale, ...
            field{level}(j).uResidual);
        far = targetSet{level}.labels == "offset 2";
        farScale = max(abs(ref.p(far)));
        fprintf(['  far-field pressure (offset 2, %d points): ', ...
            'reference scale %.3e, max |dp| %.3e, ', ...
            'far-field relative %.3e\n'], ...
            nnz(far),farScale,max(dp(far)),max(dp(far))/farScale);
        if level == 1
            for name = unique(targetSet{level}.labels).'
                if name == "offset 2"
                    continue
                end
                mask = targetSet{level}.labels == name;
                fprintf('  pressure %s: max |dp| %.3e, global relative %.3e\n', ...
                    name,max(dp(mask)),max(dp(mask))/pScale);
            end
        end
    end
end

refChangeP = max(abs(field{1}(refIndex-1).p-field{1}(refIndex).p));
refChangeT = max(vecnorm( ...
    field{1}(refIndex-1).t-field{1}(refIndex).t,2,2));
imageChangeP = max(abs(field{1}(lastImageIndex-1).p- ...
    field{1}(lastImageIndex).p));
imageChangeT = max(vecnorm(field{1}(lastImageIndex-1).t- ...
    field{1}(lastImageIndex).t,2,2));
solverChangeP = max(abs(field{1}(baseImageIndex).p- ...
    field{1}(looseIndex).p));
solverChangeT = max(vecnorm(field{1}(baseImageIndex).t- ...
    field{1}(looseIndex).t,2,2));
coarseP = max(abs(field{1}(baseImageIndex).p-field{1}(refIndex).p));
fineP = max(abs(field{2}(baseImageIndex).p-field{2}(refIndex).p));
coarseT = max(vecnorm(field{1}(baseImageIndex).t- ...
    field{1}(refIndex).t,2,2));
fineT = max(vecnorm(field{2}(baseImageIndex).t- ...
    field{2}(refIndex).t,2,2));
samplingChangeP = abs(fineP-coarseP);
samplingChangeT = abs(fineT-coarseT);
fprintf(['Refinement changes (absolute max): reference p %.3e, t %.3e; ', ...
    'image p %.3e, t %.3e; GMRES p %.3e, t %.3e; ', ...
    'sampling maxima p %.3e, t %.3e\n'], ...
    refChangeP,refChangeT,imageChangeP,imageChangeT, ...
    solverChangeP,solverChangeT,samplingChangeP,samplingChangeT);
if refChangeP > coarseP || refChangeT > coarseT
    fprintf(['Conclusion: reference-limited. The last S-only refinement ', ...
        'changes pressure or traction more than the image/reference ', ...
        'difference; agreement is consistent but not independently ', ...
        'resolved.\n']);
else
    fprintf(['Conclusion: compare the field differences above with all ', ...
        'four refinement changes before declaring resolved agreement.\n']);
end

% Display the standard paper image count against the refined S-only solve.
dp = abs(field{1}(baseImageIndex).p-field{1}(refIndex).p);
dt = vecnorm(field{1}(baseImageIndex).t- ...
    field{1}(refIndex).t,2,2);
pFloor = eps(max(1,max(dp)));
tFloor = eps(max(1,max(dt)));
fig = figure('Name','Two-sphere pressure and traction differences');
layout = tiledlayout(1,2);
ax = nexttile(layout);
points = targetSet{1}.exterior;
scatter3(ax,points(:,1),points(:,2),points(:,3), ...
    8,log10(max(dp,pFloor)),'filled');
axis(ax,'equal'); grid(ax,'on'); colorbar(ax);
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
title(ax,'log_{10}|p_{image}-p_{reference}|');
ax = nexttile(layout);
points = targetSet{1}.surface;
scatter3(ax,points(:,1),points(:,2),points(:,3), ...
    8,log10(max(dt,tFloor)),'filled');
axis(ax,'equal'); grid(ax,'on'); colorbar(ax);
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
title(ax,'log_{10}||t_{image}-t_{reference}||_2');
exportgraphics(fig,fullfile(root,'visuals', ...
    'oct17_two_sphere_pressure_traction_difference.png'),'Resolution',180);
