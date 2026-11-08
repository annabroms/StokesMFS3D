function offsets = stagnant_cap_cluster_offsets(count,extent,clustering,ratio,sigma)
%STAGNANT_CAP_CLUSTER_OFFSETS Share outer-to-inner source and surface grading.
% Geometric: extent*ratio.^k. Tapered: extent*exp(-sigma*(sqrt(N)-sqrt(N-k))).
% Here k=0:N-1; extent is a distance for sources and an angle for surface rings.
% The tapered rule follows Herremans, Huybrechs and Trefethen (2023), eq. (2.2).
% https://doi.org/10.1137/23M1551821
if nargin==0
    self_test();
    offsets = [];
    return
end
level = 0:count-1;
switch clustering
    case "geometric"
        offsets = extent*ratio.^level;
    case "tapered"
        offsets = extent*exp(-sigma*(sqrt(count)-sqrt(count-level)));
    otherwise
        error('stagnant_cap_cluster_offsets:BadClustering', ...
            'Clustering must be geometric or tapered.');
end
assert(all(isfinite(offsets)) && all(offsets>0) && all(diff(offsets)<0), ...
    'stagnant_cap_cluster_offsets:UnresolvedOffsets', ...
    'Cluster offsets must be positive and distinct; reduce clustering strength.');
end

function self_test()
% Check endpoints, ordering, tapering, and exact geometric compatibility.
geometric = stagnant_cap_cluster_offsets(5,0.35,"geometric",0.1,4);
assert(isequal(geometric,0.35*0.1.^(0:4)));
tapered = stagnant_cap_cluster_offsets(5,0.35,"tapered",0.1,4);
expected = 0.35*exp(-4*(sqrt(5)-sqrt(5:-1:1)));
assert(isequal(tapered,expected) && tapered(1)==0.35);
assert(abs(tapered(end)-0.002493434916709918)<1e-15);
assert(all(diff(diff(log(tapered)))<0));
for clustering = ["geometric","tapered"]
    assert(isempty(stagnant_cap_cluster_offsets(0,0.35,clustering,0.7,4)));
    assert(stagnant_cap_cluster_offsets(1,0.35,clustering,0.7,4)==0.35);
end
fprintf('stagnant_cap_cluster_offsets self-test passed.\n');
end
