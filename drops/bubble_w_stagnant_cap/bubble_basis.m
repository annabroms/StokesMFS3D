function [basis,coefficientCount] = bubble_basis(geometry,opt)
%BUBBLE_BASIS Build exterior basis blocks on interior source locations only.
% Columns are block ordered, then node-interleaved [x1;y1;z1;x2;...].
if nargin == 0
    bubble_self_test();
    basis = [];
    coefficientCount = [];
    return
end
baseCount = geometry.baseSourceCount;
conePoints = geometry.sourcePoints(baseCount+1:end,:);
coneCount = size(conePoints,1);
assert(coneCount==opt.n_level*opt.n_az_source);
block = struct('family',"S", ...
    'points',geometry.sourcePoints(1:baseCount,:), ...
    'normalDirections',zeros(0,3),'scale',1, ...
    'derivativeDistances',zeros(0,1));
basis = block;
if coneCount>0
    ray = conePoints-opt.centre;
    ray = ray./vecnorm(ray,2,2);
    for family = opt.refinementFamilies
        block.family = family;
        block.points = conePoints;
        block.normalDirections = zeros(0,3);
        block.derivativeDistances = zeros(0,1);
        block.scale = 1;
        if ismember(family,["T","T_t","T_n","T_n2"])
            block.normalDirections = image_stresslet_directions( ...
                ray,opt.use_orth_dir_images,family);
        elseif family=="Q_s"
            block.normalDirections = ray;
            block.scale = opt.radius^2;
        end
        basis(end+1,1) = block; %#ok<AGROW>
    end
    if opt.derivative_direction~="none"
        psi = atan2(conePoints(:,2)-opt.centre(2), ...
            conePoints(:,1)-opt.centre(1));
        if opt.derivative_direction=="radial"
            direction = [sin(opt.phi)*cos(psi), ...
                sin(opt.phi)*sin(psi), ...
                cos(opt.phi)*ones(coneCount,1)];
        else
            direction = [cos(opt.phi)*cos(psi), ...
                cos(opt.phi)*sin(psi), ...
                -sin(opt.phi)*ones(coneCount,1)];
        end
        block.family = "S_y";
        block.points = conePoints;
        block.normalDirections = direction;
        block.derivativeDistances = opt.radius- ...
            vecnorm(conePoints-opt.centre,2,2);
        block.scale = 1;
        assert(all(block.derivativeDistances>0));
        basis(end+1,1) = block;
    end
end
coefficientCount = 3*sum(arrayfun(@(b) size(b.points,1),basis));
assert(coefficientCount==3*(baseCount+coneCount* ...
    (numel(basis)-1)));
end

function bubble_self_test()
% Verify ordered cone blocks, orientations, scales, and moment contributions.
o = bubble_options(struct('sourceRequest',50,'radius',2,'phi',pi/4, ...
    'n_level',1,'n_az_source',6,'n_az_colloc',8, ...
    'refinementFamilies',"SRT_tT_nT_n2DQ_s", ...
    'derivative_direction',"cone_normal"));
g = bubble_geometry(o);
[basis,count] = bubble_basis(g,o);
assert(isequal([basis.family], ...
    ["S","S","R","T_t","T_n","T_n2","D","Q_s","S_y"]));
assert(count==3*(g.baseSourceCount+8*g.coneSourceCount));
assert(all(arrayfun(@(b) isequal(b.points,g.sourcePoints( ...
    g.baseSourceCount+1:end,:)),basis(2:end))));
assert(basis(8).scale==4);
for k = 4:6
    assert(max(abs(vecnorm(basis(k).normalDirections,2,2)-1))<1e-12);
end
assert(max(abs(sum(basis(4).normalDirections.* ...
    basis(5).normalDirections,2)))<1e-12);
assert(max(abs(basis(end).derivativeDistances- ...
    (o.radius-vecnorm(basis(end).points-o.centre,2,2))))<1e-12);
body.basis = basis;
[moments,~] = image_mobility_maps(body,o.centre);
offsets = 3*cumsum([0;arrayfun(@(b) size(b.points,1),basis(:))]);
stokesletColumn = offsets(2)+1;
rotletColumn = offsets(3)+2;
probe = zeros(count,1);
probe(stokesletColumn) = 2;
probe(rotletColumn) = 3;
assert(norm(moments(1:3,:)*probe-[2;0;0])<1e-13);
expectedTorque = cross(basis(2).points(1,:)-o.centre,[2 0 0]).' ...
    +[0;3;0];
assert(norm(moments(4:6,:)*probe-expectedTorque)<1e-12);
fprintf('bubble_basis self-test passed.\n');
end
