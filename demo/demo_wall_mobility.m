%DEMO_WALL_MOBILITY Mobility and resistance round trip for five equal-height
% spheres above a no-slip plane. Random force/torque data drive the mobility
% solve; its rigid velocities are then prescribed in resistance. The script
% reports the independent surface residuals and dimensionless force round-trip
% error using direct Blake evaluation and the wall-aware one-body preconditioner.

setpath;
P = 5;
a = 1;
mu = 1;
[q,~] = grow_cluster_wall(P,a,1,0.5,20261026,"equalHeight");
[rin,rout,opt] = init_spheres_wall(q,a,struct( ...
    'N',500,'mu',mu,'backend',"blake",'preconditioner',"wall"));
opt.preparedProblem = prepare_wall_problem(q,rin,rout,opt);

FT = randn(6*P,1); 
[U,~,~,mobilitySurfaceError] = ...
    solve_mobility_wall(q,rin,rout,FT,opt);
[FTback,~,~,resistanceSurfaceError] = ...
    solve_resistance_wall(q,rin,rout,U,opt);

forceUnits = repmat([mu*a*ones(3,1);mu*a^2*ones(3,1)],P,1);
twoWayError = norm((FTback-FT)./forceUnits,inf) / ...
	max(norm(FT./forceUnits,inf),realmin);
fprintf('Mobility surface residual:   %.3e\n',mobilitySurfaceError);
fprintf('Resistance surface residual: %.3e\n',resistanceSurfaceError);
fprintf('Force round-trip error:      %.3e\n',twoWayError);