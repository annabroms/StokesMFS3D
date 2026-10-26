function validation = bubble_validation_grid(opt,maxCollocation,finestOffset)
%BUBBLE_VALIDATION_GRID Independent shifted surface grid and edge rings.
if nargin == 0
    bubble_self_test();
    validation = [];
    return
end
assert(isfield(opt,'geometryOptions'),'Complete options with bubble_options.');
validation = stagnant_cap_validation_grid( ...
    opt.geometryOptions,maxCollocation,finestOffset);
assert(validation.count>=5*maxCollocation);
if opt.phi>0 && opt.phi<pi
    assert(any(validation.cap) && any(validation.mobile));
    assert(validation.minimumEdgeDistance>0);
end
end

function bubble_self_test()
% Check independent density, both edge sides, and endpoint handling.
o = bubble_options(struct('sourceRequest',50,'ringLevels',2));
g = bubble_geometry(o);
v = bubble_validation_grid(o,g.collocationCount,g.ringOffsets(end));
assert(v.count>=5*g.collocationCount);
assert(any(v.cap & v.nearEdge) && any(v.mobile & v.nearEdge));
assert(v.minimumEdgeDistance<g.minimumEdgeDistance);
o = bubble_options(struct('sourceRequest',50,'phi',0));
g = bubble_geometry(o);
v = bubble_validation_grid(o,g.collocationCount,0);
assert(~any(v.cap) && v.count>=5*g.collocationCount);
fprintf('bubble_validation_grid self-test passed.\n');
end
