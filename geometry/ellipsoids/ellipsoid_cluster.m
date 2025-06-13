function [E R t xnear] = ellipsoid_cluster(E0,K,dmin)
% ELLIPSOID_CLUSTER   grow random dmin-close cluster ellipsoids, same semiaxes
%
% [E R t xnear] = ellipsoid_cluster(E0,K,dmin) returns cell arrays defining K
%  ellipsoids of random rotations and locations, but achieving close to dmin
%  distance between each an at least one other.
%  E0 = [a b c] semiaxes of all ellipsoids (same shape).
%  Output jth ellipsoid defined by semiaxes E{j}, rotation 3x3 matrix R{j},
%  and translation center t{j}. xnear{j} lists point of nearest touching
%  on the jth ellipsoid, as a 3x1 vector (apart from j=1 which is empty).
%
% Without arguments, does self test with pic

% Barnett 8/15/21
if nargin==0, test_ellipsoid_cluster; return; end

% Euler rot mats
Rz = @(t) [cos(t) -sin(t) 0;sin(t) cos(t) 0; 0 0 1];  % z-ax rot mat
Ry = @(t) [cos(t) 0 -sin(t); 0 1 0; sin(t) 0 cos(t)];  % y-ax rot mat

t{1} = zeros(3,1); 
R{1} = eye(3); 

%Try with a random direction
u = randn(3,1); u=u/norm(u);
alpha = atan2(u(2),u(1));
beta = acos(u(3));
gamma = 2*pi*rand;
R{1} = Rz(alpha) * Ry(beta) * Rz(gamma);

E{1} = E0; xnear{1} = [];   % object 1
for j=2:K
  fprintf('starting ellipsoid j=%d...\n',j)
  fail = 1; 
  while fail
  u = randn(3,1); u=u/norm(u);  % rand on S2 defines (alpha, beta) Euler
  alpha = atan2(u(2),u(1));
  beta = acos(u(3));
  gamma = 2*pi*rand;
  R{j} = Rz(alpha) * Ry(beta) * Rz(gamma);  % Euler rot (unif on SO(3) I think)
  %R{j} = eye(3); %debug!
  %R{j} = R{1};
  E{j} = E0;           % all same shape
  dir = randn(3,1); dir=dir/norm(dir);  % direction to approach from, rand on S2
  distprev = @(s) min_dist_ellipsoids(E{j},R{j},s*dir,E(1:j-1),R(1:j-1),t(1:j-1));
  % old diagnosis code
  %figure; s=2:0.1:10; for i=1:numel(s), dd(i)=distprev(s(i)); end % diagnose
  %plot(s,dd, '+-'); title(sprintf('j=%d',j)); drawnow
  %s = fzero(@(s) distprev(s)-dmin,[1 10]);   % try to hit dmin - fails, too close
  
  dsteps = 10.^[1:-1:-4];   % careful approach so distprev not eval too close!
  s = 2+2*max(E0);   % start dist far enough. we assume semiaxes O(1)
  s = 2+6*max(E0);
  s = 2+12*max(E0);
  for dstep = dsteps
    d=distprev(s)-dmin;      % aiming this to hit 0 from above
    while d>dstep
      s=s-dstep; d=distprev(s)-dmin; %[s,d]     % chug down
    end
  end
  [d,x,y,i] = min_dist_ellipsoids(E{j},R{j},s*dir,E(1:j-1),R(1:j-1),t(1:j-1));
  if d<1e-3
    %  disp('redo, too close')
  elseif abs(d-dmin)/dmin<1e-2
%       disp('Final d')
%       d
      fail = 0;
  else
  
      %disp('redo');
  end
end
  fprintf('\tellipsoid j=%d achieved dist %g from i=%d\n',j,d,i)
  t{j} = s*dir;
  xnear{j} = x;
end


%%%%%%%%%%%%%
function test_ellipsoid_cluster
E0 = [.5,1,1.5];     % baseline semiaxes
E0 = [.5 .5 1];
K=10;   % build cluster of K same shape ellipsoids
dmin = 0.01;    % target dist of each to prev ellipsoids
rng(0);    % seed
tic; [E R t xnear] = ellipsoid_cluster(E0,K,dmin); toc     % do it

figure; for j=1:K                 % show it
  colorvec = mod(j*[.42,.1,.29],1);   % some irrationals so colors unique
  %Nvec = [30,60]; % I don't know what to put here... /Anna 
  show_ellipsoid_v2(E{j},R{j},t{j},colorvec); hold on;
  x = xnear{j}; if ~isempty(x), plot3(x(1),x(2),x(3),'k.','markersize',20); end
  text(t{j}(1),t{j}(2),t{j}(3),sprintf('%d',j), 'fontsize',20);
end