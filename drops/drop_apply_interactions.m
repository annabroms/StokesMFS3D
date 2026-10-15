function interaction = drop_apply_interactions(lambda, innerSources, ...
    surfacePoints, normals, factors, opt)
%DROP_APPLY_INTERACTIONS Apply all off-diagonal recompleted drop blocks.
% Inputs lambda{j} are force-free, node-interleaved exterior strengths.
% The output has body-major [continuity;normal velocity;tangent traction]
% rows, already scaled exactly like the local one-body matrices.
% Direct source-target batches can later be replaced by an accelerated sum.

if nargin == 0
    self_test();
    return
end

bodyCount = numel(innerSources);
assert(numel(lambda)==bodyCount && numel(surfacePoints)==bodyCount && ...
    numel(normals)==bodyCount && numel(factors)==bodyCount);
targetCounts = cellfun(@(points) size(points,1),surfacePoints);
interaction = zeros(6*sum(targetCounts),1);
rowOffset = 0;

for targetBody = 1:bodyCount
    targetCount = targetCounts(targetBody);
    targetPoints = surfacePoints{targetBody};
    targetNormals = normals{targetBody};
    tangentOne = factors{targetBody}.tangentOne;
    tangentTwo = factors{targetBody}.tangentTwo;
    assert(isequal(size(targetNormals),size(targetPoints)) && ...
        isequal(size(tangentOne),size(targetPoints)) && ...
        isequal(size(tangentTwo),size(targetPoints)));
    normalVelocity = zeros(targetCount,1);
    tangentTraction = zeros(targetCount,2);

    for sourceBody = 1:bodyCount
        if sourceBody == targetBody
            continue
        end
        sources = innerSources{sourceBody};
        sourceCount = size(sources,1);
        assert(numel(lambda{sourceBody}) == 3*sourceCount);
        strengths = reshape(lambda{sourceBody},3,sourceCount).';

        for firstTarget = 1:opt.interactionBatchSize:targetCount
            rows = firstTarget:min(firstTarget+opt.interactionBatchSize-1, ...
                targetCount);
            targets = targetPoints(rows,:);
            targetNormalBatch = targetNormals(rows,:);
            rx = targets(:,1)-sources(:,1).';
            ry = targets(:,2)-sources(:,2).';
            rz = targets(:,3)-sources(:,3).';
            distanceSquared = rx.^2+ry.^2+rz.^2;
            assert(all(distanceSquared > 0,'all'), ...
                'A cross-drop source coincides with a target.');
            inverseDistance = 1./sqrt(distanceSquared);
            radialStrength = rx.*strengths(:,1).' + ...
                ry.*strengths(:,2).' + rz.*strengths(:,3).';
            velocityRadial = radialStrength.*inverseDistance.^3;
            velocity = [ ...
                sum(strengths(:,1)'.*inverseDistance+rx.*velocityRadial,2), ...
                sum(strengths(:,2)'.*inverseDistance+ry.*velocityRadial,2), ...
                sum(strengths(:,3)'.*inverseDistance+rz.*velocityRadial,2)] ...
                /(8*pi);
            normalVelocity(rows) = normalVelocity(rows) + ...
                sum(velocity.*targetNormalBatch,2);

            radialNormal = rx.*targetNormalBatch(:,1) + ...
                ry.*targetNormalBatch(:,2) + ...
                rz.*targetNormalBatch(:,3);
            % The supplied drop formulation has -3/(4*pi), opposite the
            % sign of the repository's getTractionDirectSLP routine.
            tractionWeight = -3/(4*pi)*radialNormal.*radialStrength.* ...
                inverseDistance.^5;
            traction = [sum(tractionWeight.*rx,2), ...
                sum(tractionWeight.*ry,2), ...
                sum(tractionWeight.*rz,2)];
            tangentTraction(rows,1) = tangentTraction(rows,1) + ...
                sum(traction.*tangentOne(rows,:),2);
            tangentTraction(rows,2) = tangentTraction(rows,2) + ...
                sum(traction.*tangentTwo(rows,:),2);
        end
    end

    interaction(rowOffset+3*targetCount+(1:targetCount)) = ...
        normalVelocity/opt.velocityScale;
    tangentTraction = (1-opt.interiorViscosity(targetBody))* ...
        opt.tractionRowWeight*tangentTraction/opt.tractionScale;
    interaction(rowOffset+4*targetCount+(1:2*targetCount)) = ...
        reshape(tangentTraction.',[],1);
    rowOffset = rowOffset+6*targetCount;
end
end

function self_test()
% Compare direct pair sums with the existing explicit interaction blocks.
axisNodes = [eye(3);-eye(3)];
corners = [1 1 1;1 1 -1;1 -1 1;-1 1 1; ...
           1 -1 -1;-1 1 -1;-1 -1 1;-1 -1 -1]/sqrt(3);
centres = [-1.6 0 0;1.6 0 0];
inner = cell(1,2); outer = cell(1,2);
surface = cell(1,2); normals = cell(1,2); factors = cell(1,2);
for body = 1:2
    inner{body} = centres(body,:)+0.6*axisNodes;
    outer{body} = centres(body,:)+2*axisNodes;
    surface{body} = centres(body,:)+corners;
    normals{body} = corners;
    factors{body} = oneBodyPrecondDrop(inner{body},outer{body}, ...
        surface{body},normals{body},2,1/(8*pi),3,5,1.7,1e-12,false);
end
lambda = {sin((1:18).'),cos((1:18).')};
for body = 1:2
    nodeValues = reshape(lambda{body},3,[]);
    lambda{body} = reshape(nodeValues-mean(nodeValues,2),[],1);
end
opt.interiorViscosity = [2 2];
opt.velocityScale = 3;
opt.tractionScale = 5;
opt.tractionRowWeight = 1.7;
opt.interactionBatchSize = 3;
actual = drop_apply_interactions(lambda,inner,surface,normals,factors,opt);
expected = zeros(96,1);
for body = 1:2
    other = 3-body;
    block = drop_interaction_block(inner{other},surface{body}, ...
        normals{body},factors{body}.tangentOne, ...
        factors{body}.tangentTwo,2);
    local = block*[lambda{other};zeros(18,1)];
    local(1:32) = local(1:32)/3;
    local(33:end) = 1.7*local(33:end)/5;
    expected((body-1)*48+(1:48)) = local;
end
assert(norm(actual-expected,inf) < 1e-12);
assert(~any(actual(1:24)) && ~any(actual(49:72)));
fprintf('drop_apply_interactions self-test passed.\n');
end
