function paths = get_streamline(startXZ, velocityAtXZ, distanceToStop, ...
    maxArcLength, maxStep, minimumSpeed, odeTolerances)
%GET_STREAMLINE Trace planar streamlines with MATLAB's adaptive ode45.
%
% INPUTS
%   startXZ         K-by-2 starting positions [x,z], one per streamline.
%   velocityAtXZ    Function returning [u_x,u_z] at a 1-by-2 position.
%   distanceToStop  Function returning a positive distance while tracing
%                   is allowed and zero at the stopping boundary. For the
%                   global plot this is the distance to the plot box; for
%                   an interior plot it is a-r from the drop surface.
%   maxArcLength    Longest distance traced from a seed in each direction.
%   maxStep         Largest arc-length step permitted to ode45.
%   minimumSpeed    Stop when the physical speed falls below this value.
%   odeTolerances   Optional [relative, absolute] ode45 tolerances; when
%                   omitted, use [1e-6, 1e-8*maxStep].
%
% OUTPUT
%   paths           K-by-1 cell array. Each cell contains an L-by-2 curve
%                   [x,z], ordered from backward end through the start to
%                   forward end. A stagnant or invalid seed stays a point.
%
% METHOD
%   ode45 solves dx/ds = +/-u(x)/|u(x)|. The variable s is distance along
%   the curve, not physical time. Normalization changes only the rate of
%   tracing, not the streamline geometry. The two directions are joined.
%
% STOPPING
%   Terminal ode45 events locate the domain boundary, weak flow, or a
%   nearly closed loop. Seeds must start strictly inside the allowed area.
%   The global distance function mentions no sphere, so its lines can
%   cross interfaces when velocityAtXZ switches physical regions.

if nargin == 0
    self_test();
    return
end
if nargin < 7
    odeTolerances = [1e-6,1e-8*maxStep];
end

assert(size(startXZ,2)==2 && numel(odeTolerances)==2 && ...
    maxArcLength>0 && maxStep>0 && minimumSpeed>=0);
paths = cell(size(startXZ,1),1);
for seedIndex = 1:size(startXZ,1)
    seed = startXZ(seedIndex,:);
    if distanceToStop(seed)<=0 || norm(velocityAtXZ(seed))<=minimumSpeed
        paths{seedIndex} = seed;
        continue
    end

    % Event location keeps the end points on the requested boundary.
    odeOptions = odeset('RelTol',odeTolerances(1), ...
        'AbsTol',odeTolerances(2),'MaxStep',maxStep, ...
        'Events',@(s,x) stopEvents(s,x,seed,velocityAtXZ, ...
        distanceToStop,minimumSpeed,maxStep));
    [~,forward,~,~,forwardEvents] = ode45( ...
        @(s,x) unitVelocity(s,x,+1,velocityAtXZ,minimumSpeed), ...
        [0 maxArcLength],seed.',odeOptions);
    if any(forwardEvents==3)
        % One complete closed loop is enough; the reverse trace duplicates it.
        paths{seedIndex} = [forward;seed];
        continue
    end
    [~,backward] = ode45( ...
        @(s,x) unitVelocity(s,x,-1,velocityAtXZ,minimumSpeed), ...
        [0 maxArcLength],seed.',odeOptions);
    paths{seedIndex} = [flipud(backward(2:end,:));forward];
end
end

function derivative = unitVelocity(~, position, signOfTrace, velocityAtXZ, minimumSpeed)
% Supply ode45 with a unit tangent to the physical velocity field.
velocity = velocityAtXZ(position.');
speed = norm(velocity);
if speed <= minimumSpeed
    derivative = zeros(2,1);
else
    derivative = signOfTrace*velocity(:)/speed;
end
end

function [value,isTerminal,direction] = stopEvents( ...
    arcLength,position,seed,velocityAtXZ,distanceToStop,minimumSpeed,maxStep)
% Stop at the allowed-region edge, stagnation, or return to the seed.
point = position.';
value = [distanceToStop(point); ...
    norm(velocityAtXZ(point))-minimumSpeed;1];
if arcLength > 10*maxStep
    value(3) = norm(point-seed)-0.6*maxStep;
end
isTerminal = [1;1;1];
direction = [-1;-1;-1];
end

function self_test()
% Straight, circular, crossing, and stagnant fields exercise ode45 events.
box = @(point) min([point(1)+1,1-point(1), ...
    point(2)+1,1-point(2)]);
straight = get_streamline([0 0],@(point) [1 0],box,3,0.05,1e-8);
assert(iscell(straight) && isscalar(straight) && ...
    size(straight{1},2)==2 && ...
    max(abs(straight{1}(:,2)))<1e-10 && ...
    min(straight{1}(:,1))<-0.99 && max(straight{1}(:,1))>0.99);

circle = get_streamline([1 0], ...
    @(point) [-point(2) point(1)],@(point) 2-norm(point), ...
    8,0.04,1e-8);
assert(size(circle{1},1)>30 && ...
    max(abs(vecnorm(circle{1},2,2)-1))<2e-3);

crossingField = @(point) [1,0.2*(point(1)>=0)];
crossing = get_streamline([-0.5 0],crossingField,box,3,0.03,1e-8);
assert(min(crossing{1}(:,1))<-0.9 && ...
    max(crossing{1}(:,1))>0.9 && max(crossing{1}(:,2))>0.05);
stagnant = get_streamline([0 0],@(point) [0 0],box,1,0.1,1e-6);
assert(isequal(stagnant{1},[0 0]));
fprintf('get_streamline self-test passed.\n');
end
