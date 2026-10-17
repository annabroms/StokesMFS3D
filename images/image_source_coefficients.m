function coeff = image_source_coefficients(basis,strength)
%IMAGE_SOURCE_COEFFICIENTS Split and scale interleaved image strengths.
%   image_source_coefficients() checks the complete block ordering.

if nargin == 0
    self_test();
    coeff = {};
    return
end
coeff = cell(numel(basis),1);
offset = 0;
for b = 1:numel(basis)
    count = size(basis(b).points,1);
    cols = offset+(1:3*count);
    coeff{b} = basis(b).scale*reshape(strength(cols),3,[]).';
    offset = offset+3*count;
end
assert(offset == numel(strength));
end

function self_test()
%SELF_TEST checks two differently scaled source blocks.
basis = repmat(struct('points',zeros(0,3),'scale',1),2,1);
basis(1).points = zeros(2,3);
basis(2).points = zeros(1,3);
basis(2).scale = 2;
coeff = image_source_coefficients(basis,(1:9).');
assert(isequal(coeff{1},[1,2,3;4,5,6]));
assert(isequal(coeff{2},[14,16,18]));
fprintf('image_source_coefficients self-test passed.\n');
end
