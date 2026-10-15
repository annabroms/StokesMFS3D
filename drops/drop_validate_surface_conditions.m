function [maximumErrors, pointErrors] = drop_validate_surface_conditions( ...
    lambda, gamma, translations, innerSources, outerSources, ...
    validationPoints, normals, activeDirections, S1, viscosity, batchSize)
%DROP_VALIDATE_SURFACE_CONDITIONS Check all three conditions off collocation.
% Output rows correspond to bodies; columns are physical maxima for
% velocity continuity, normal velocity, and tangential traction. Optional
% pointErrors{body} has the same three physical residuals at every
% validation point, in its input order, for surface plots.

if nargin == 0
    self_test();
    return
end

bodyCount = numel(innerSources);
assert(numel(lambda)==bodyCount && numel(gamma)==bodyCount && ...
    numel(outerSources)==bodyCount && numel(validationPoints)==bodyCount && ...
    numel(normals)==bodyCount && ...
    isequal(size(translations),[bodyCount 3]) && ...
    isequal(size(activeDirections),[bodyCount 3]));
if isscalar(S1), S1 = repmat(S1,bodyCount,1); end
if isscalar(viscosity), viscosity = repmat(viscosity,bodyCount,1); end
assert(numel(S1)==bodyCount && numel(viscosity)==bodyCount && ...
    all(viscosity>0) && batchSize>=1 && batchSize==round(batchSize));
maximumErrors = zeros(bodyCount,3);
if nargout > 1
    pointErrors = cell(1,bodyCount);
end

for body = 1:bodyCount
    pointCount = size(validationPoints{body},1);
    assert(isequal(size(normals{body}),[pointCount 3]));
    if nargout > 1
        pointErrors{body} = zeros(pointCount,3);
    end
    direction = activeDirections(body,:).';
    direction = direction/norm(direction);
    for firstPoint = 1:batchSize:pointCount
        selected = firstPoint:min(firstPoint+batchSize-1,pointCount);
        targets = validationPoints{body}(selected,:);
        targetNormals = normals{body}(selected,:);

        selfVelocity = drop_velocity_block( ...
            innerSources{body},targets,1)*lambda{body};
        interiorVelocity = drop_velocity_block( ...
            outerSources{body},targets,viscosity(body))*gamma{body};
        incidentVelocity = zeros(size(selfVelocity));
        for other = 1:bodyCount
            if other~=body
                incidentVelocity = incidentVelocity + ...
                    drop_velocity_block(innerSources{other},targets,1)* ...
                    lambda{other};
            end
        end
        continuityError = vecnorm(reshape( ...
            selfVelocity-interiorVelocity,3,[]).',2,2);
        exteriorVelocity = reshape(selfVelocity+incidentVelocity,3,[]).';
        normalError = abs(sum(exteriorVelocity.*targetNormals,2)- ...
            targetNormals*translations(body,:).');

        tractionJump = drop_traction_block( ...
            innerSources{body},targets,targetNormals)*lambda{body} - ...
            drop_traction_block( ...
            outerSources{body},targets,targetNormals)*gamma{body};
        for other = 1:bodyCount
            if other~=body
                tractionJump = tractionJump + ...
                    (1-viscosity(body))*drop_traction_block( ...
                    innerSources{other},targets,targetNormals)*lambda{other};
            end
        end
        tractionJump = reshape(tractionJump,3,[]).';
        tangentJump = tractionJump- ...
            sum(tractionJump.*targetNormals,2).*targetNormals;
        normalDotDirection = targetNormals*direction;
        prescribedTraction = S1(body)*normalDotDirection.* ...
            (direction.'-normalDotDirection.*targetNormals);
        tractionError = vecnorm(tangentJump-prescribedTraction,2,2);
        batchErrors = [continuityError,normalError,tractionError];
        maximumErrors(body,:) = max(maximumErrors(body,:), ...
            max(batchErrors,[],1));
        if nargout > 1
            pointErrors{body}(selected,:) = batchErrors;
        end
    end
end
end

function self_test()
% Zero data satisfy every condition; a radial Stokeslet checks magnitudes.
axisNodes = [eye(3);-eye(3)];
normals = {[1 0 0;0 1 0;0 0 1]};
points = normals;
lambda = {zeros(18,1)};
gamma = {zeros(18,1)};
[errors,pointErrors] = drop_validate_surface_conditions( ...
    lambda,gamma,zeros(1,3), ...
    {0.6*axisNodes},{2*axisNodes},points,normals,[0 0 1],0,2,2);
assert(isequal(size(errors),[1 3]) && ~any(errors));
assert(isequal(size(pointErrors{1}),[3 3]) && ~any(pointErrors{1},'all'));
lambda{1}(1) = 1;
[errors,pointErrors] = drop_validate_surface_conditions( ...
    lambda,gamma,zeros(1,3), ...
    {0.6*axisNodes},{2*axisNodes},points,normals,[0 0 1],0,2,2);
assert(all(isfinite(errors)) && all(errors>0));
assert(norm(errors-max(pointErrors{1},[],1),inf)<1e-14);
fprintf('drop_validate_surface_conditions self-test passed.\n');
end
