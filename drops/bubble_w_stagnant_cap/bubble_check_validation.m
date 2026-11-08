function check = bubble_check_validation( ...
    solution,opt,validation,fineValidation)
%BUBBLE_CHECK_VALIDATION Compare independent grids at two edge resolutions.
if nargin == 0
    bubble_self_test();
    check = [];
    return
end
geometry = solution.grid;
if opt.phi>0 && opt.phi<pi
    if isempty(geometry.ringOffsets)
        finestOffset = opt.edgeBandAngle/2;
    else
        finestOffset = geometry.ringOffsets(end);
    end
else
    finestOffset = 0;
end
if nargin<3
    validation = bubble_validation_grid( ...
        opt,geometry.collocationCount,finestOffset);
else
    assert(nargin==4, ...
        'Supply both shared validation grids or neither.');
end
assert(validation.count>=5*geometry.collocationCount);
check.validation = validation;
check.report = bubble_validate(solution,check.validation,opt);
fineOffset = 0;
if opt.phi>0 && opt.phi<pi
    fineOffset = min(finestOffset/2, ...
        check.validation.minimumEdgeDistance/opt.radius);
end
if nargin<3
    fineOptions = opt;
    fineOptions.validationBaseRequest = max(opt.validationBaseRequest, ...
        2*geometry.collocationCount);
    fineOptions = bubble_options(fineOptions);
    fineValidation = bubble_validation_grid( ...
        fineOptions,2*geometry.collocationCount,fineOffset);
end
assert(fineValidation.count>validation.count);
if opt.phi>0 && opt.phi<pi
    assert(fineValidation.minimumEdgeDistance< ...
        validation.minimumEdgeDistance);
end
check.fineValidation = fineValidation;
check.fineReport = bubble_validate(solution,check.fineValidation,opt);
first = check.report.maxima;
second = check.fineReport.maxima;
denominator = max(max(abs(first),abs(second)),1e-12);
check.relativeChange = abs(second-first)./denominator;
check.unresolved = any(check.relativeChange(isfinite( ...
    check.relativeChange))>0.1);
if opt.phi>0 && opt.phi<pi
    fprintf(['Validation %d points, closest edge %.4g a; refined %d ', ...
        'points, closest edge %.4g a.\n'],check.validation.count, ...
        check.validation.minimumEdgeDistance/opt.radius, ...
        check.fineValidation.count, ...
        check.fineValidation.minimumEdgeDistance/opt.radius);
else
    fprintf('Endpoint validation: %d points; refined %d points; no cap edge.\n', ...
        check.validation.count,check.fineValidation.count);
end
for k = 1:4
    fprintf('%s: %.4g -> %.4g (relative change %.3g).\n', ...
        check.report.conditionNames(k),first(k),second(k), ...
        check.relativeChange(k));
end
if check.unresolved
    fprintf('Validation maxima changed by more than 10%%; edge resolution remains unresolved.\n');
end
end

function bubble_self_test()
% Check the finer grid and inapplicable endpoint conditions.
o = bubble_options(struct('sourceRequest',50,'phi',pi/3, ...
    'ringLevels',2,'batchSize',128));
s = bubble_solve(o,[0;0;-1]);
c = bubble_check_validation(s,o);
assert(c.validation.count>=5*s.grid.collocationCount);
assert(c.fineValidation.count>c.validation.count);
assert(c.fineValidation.minimumEdgeDistance< ...
    c.validation.minimumEdgeDistance);
assert(all(isfinite(c.relativeChange)));
shared = bubble_check_validation(s,o,c.validation,c.fineValidation);
assert(isequal(shared.report.maxima,c.report.maxima) && ...
    isequal(shared.relativeChange,c.relativeChange));
o = bubble_options(struct('sourceRequest',50,'phi',0));
s = bubble_solve(o,[0;0;-1]);
c = bubble_check_validation(s,o);
assert(isnan(c.report.maxima(1)) && isnan(c.report.maxima(3)));
fprintf('bubble_check_validation self-test passed.\n');
end
