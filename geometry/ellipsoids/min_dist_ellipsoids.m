function [d x y i] = min_dist_ellipsoids(E,R,t,Es,Rs,ts)
% return min dist of one ellipsoid (E,R,t) to set of other ellipsoids
%
% [d x y i] = min_dist_ellipsoids(E,R,t,Es,Rs,ts)
% E = semiaxes list, R = 3x3 rot mat, t = translation center.
%  Es, Rs, ts are cell arrays of same for other ellipsoids.
% Outputs:
% d min dist
% x,y are pair of points achieving the min dist
% i index of ellipsoid in list that triggered it

% Barnett 8/15/21
d = inf;
d = 1e6; 
for i=1:numel(Es)
  [di xi yi] = dist_ellipsoids(E,R,t,Es{i},Rs{i},ts{i});
  if di<d
    x = xi; y = yi; d = di;
  end
end