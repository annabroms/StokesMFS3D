function test_bubble_clustering()
%TEST_BUBBLE_CLUSTERING Verify independent grading without running experiments.
stagnant_cap_cluster_offsets();
defaults = bubble_options(struct());
assert(defaults.sourceClustering=="geometric" && ...
    defaults.collocClustering=="geometric");
input = struct('sourceRequest',12,'phi',pi/4,'n_level',3, ...
    'ringLevels',6,'n_az_source',8,'n_az_colloc',12, ...
    'innerSourceStartDistance',0.35,'sourceClusterRatio',0.7, ...
    'collocAngularExtent',0.2,'collocClusterRatio',0.5);
baseOpt = bubble_options(input);
base = bubble_geometry(baseOpt);
assert(isequal(base.ringOffsets,0.2*0.5.^(0:5).'));
for sourceMode = ["geometric","tapered"]
    for collocMode = ["geometric","tapered"]
        input.sourceClustering = sourceMode;
        input.collocClustering = collocMode;
        o = bubble_options(input);
        assert(isequal(o,bubble_options(o)));
        g = bubble_geometry(o);
        assert(g.sourceCount==base.sourceCount && ...
            g.collocationCount==base.collocationCount);
        assert(isequal(g.sourcePoints(1:g.baseSourceCount,:), ...
            base.sourcePoints(1:base.baseSourceCount,:)));
        if sourceMode=="geometric"
            depths = 0.35*0.7.^(0:2);
            assert(isequal(g.sourcePoints,base.sourcePoints));
        else
            depths = 0.35*exp(-4*(sqrt(3)-sqrt(3:-1:1)));
        end
        actual = 1-vecnorm(g.sourcePoints(g.baseSourceCount+1:end,:),2,2);
        assert(max(abs(actual-repelem(depths(:),8)))<1e-13);
        if collocMode=="geometric"
            offsets = 0.2*0.5.^(0:5);
            assert(isequal(g.points,base.points));
        else
            offsets = 0.2*exp(-4*(sqrt(6)-sqrt(6:-1:1)));
        end
        assert(isequal(g.ringOffsets,offsets(:)));
        for level = 1:6
            rows = g.collocationRingIndex==level;
            assert(max(abs(abs(g.theta(rows)-o.phi)-offsets(level)))<1e-13);
            assert(nnz(rows & g.cap)==g.ringAzimuthCount && ...
                nnz(rows & g.mobile)==g.ringAzimuthCount);
        end
    end
end
% The shared two-fluid builder must grade both source sides consistently.
raw = stagnant_cap_refined_geometry(o.geometryOptions);
assert(max(abs(raw.sourceDistances.inner-depths(:)))<1e-14);
outerDepths = o.geometryOptions.outerSourceStartDistance*depths(:)/0.35;
assert(max(abs(raw.sourceDistances.outer-outerDepths))<1e-14);
% Tapering surface rings must also work without enhancing sources.
uniform = bubble_options(struct('sourceRequest',12,'n_level',0, ...
    'ringLevels',4,'collocClustering',"tapered",'collocClusterSigma',3));
g0 = bubble_geometry(uniform);
expected = uniform.edgeBandAngle*exp(-3*(sqrt(4)-sqrt(4:-1:1)));
assert(isequal(g0.ringOffsets,expected(:)) && g0.coneSourceCount==0);
for phi = [0 pi]
    endpoint = uniform;
    endpoint.phi = phi;
    endpoint.edgeBandAngle = 0;
    gEnd = bubble_geometry(bubble_options(endpoint));
    assert(isempty(gEnd.ringOffsets));
end
for name = ["sourceClustering","collocClustering"]
    try
        bubble_options(struct(name,"invalid"));
        error('test_bubble_clustering:AcceptedMode','Invalid mode was accepted.');
    catch ME
        assert(startsWith(ME.identifier,'stagnant_cap_options:Bad'));
    end
end
% Exercise the new geometry with both previously added solver scalings.
o.useColumnScaling = true;
o.useAreaWeights = true;
s = bubble_solve(o,[0;0;-1]);
assert(s.forceResidual<1e-9 && all(isfinite(s.lambda)));
assert(abs(sum(s.collocationAreas)-4*pi)<1e-10);
fprintf('test_bubble_clustering passed: all grading pairs, endpoints, and scaled solve.\n');
end
