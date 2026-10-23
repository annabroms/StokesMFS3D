function test_drop_image_source_directions()
%TEST_DROP_IMAGE_SOURCE_DIRECTIONS Check three active-direction configurations.
% The selected source axes depend on image-ray geometry, not active e.
q = [-1.05 0 0;1.05 0 0]; a = ones(2,1); S1 = 2*ones(2,1);
opt = drop_image_options(struct('N',30,'collocationRatio',2, ...
    'N_image',1,'imageFamilies',"SDRT_tT_nT_n2Q_s", ...
    'mirror_ext',true,'capPointCounts',[12 12], ...
    'interiorViscosity',10,'svdCutoff',1e-12, ...
    'gmresTolerance',1e-10,'maxIterations',100),2,a);
directions = cat(3,[1 0 0;1 0 0],[0 0 1;0 0 1], ...
    [1 2 3;2 -1 1]);
translations = zeros(2,3,3);
for k = 1:3
    e = directions(:,:,k);
    [U,lambda,gamma,g,info] = solve_drops_with_images(q,a,S1,e,opt);
    [x,n] = drop_image_validation_grid(g,50,[21 11],[pi/5 pi/60]);
    error = drop_image_validate_surface_conditions( ...
        g,lambda,gamma,U,S1,e,info.options,x,n);
    assert(all(isfinite(U),'all') && all(isfinite(error),'all'));
    assert(isequal(image_family_list(opt.imageFamilies), ...
        ["S","D","R","T_t","T_n","T_n2","Q_s"]));
    ray = g.body(1).imageDirections;
    basis = g.body(1).basis;
    assert(norm(basis(5).normalDirections-ray,'fro')<1e-14);
    assert(max(abs(sum(basis(6).normalDirections.*ray,2)))<1e-14);
    assert(max(abs(sum(basis(7).normalDirections.*ray,2)))<1e-14);
    assert(max(abs(sum(basis(6).normalDirections.* ...
        basis(7).normalDirections,2)))<1e-14);
    translations(:,:,k) = U;
end
assert(norm(translations(:,:,1)-translations(:,:,2),'fro')>1e-5);
assert(norm(translations(:,:,2)-translations(:,:,3),'fro')>1e-5);
fprintf('test_drop_image_source_directions passed for axial, transverse, and distinct oblique e.\n');
end
