function [directions,weights] = stagnant_cap_design(minimumCount)
%STAGNANT_CAP_DESIGN Return the first available spherical design of at least N.
if nargin == 0
    self_test();
    return
end
validateattributes(minimumCount,{'numeric'}, ...
    {'real','finite','scalar','integer','positive'});
catalogue = fullfile(fileparts(which('get_sphdesign')), ...
    'spherical_designs','filelist.txt');
fid = fopen(catalogue,'r');
assert(fid>=0,'Spherical-design catalogue is unavailable.');
data = textscan(fid,'%*[^.].%d');
fclose(fid);
counts = data{1};
index = find(counts>=minimumCount,1);
if isempty(index)
    error('stagnant_cap_design:CatalogueExceeded', ...
        'Requested %d points; largest spherical design has %d.', ...
        minimumCount,max(counts));
end
[directions,weights] = get_sphdesign(counts(index));
assert(size(directions,1)>=minimumCount);
end

function self_test()
%SELF_TEST Check catalogue selection, sphere radius, and quadrature mass.
[x,w] = stagnant_cap_design(51);
assert(size(x,1)>=51 && size(x,2)==3);
assert(max(abs(vecnorm(x,2,2)-1))<1e-12);
assert(abs(sum(w)-4*pi)<1e-12);
fprintf('stagnant_cap_design self-test passed.\n');
end
