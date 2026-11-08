function families = ring_family_list(value)
%RING_FAMILY_LIST Parse ring family combinations, including S_y.
% Input: character vector, string scalar or string array of family tokens.
% Output: row string array in supplied order; empty input is the GL baseline.
% S_y is the depth-scaled source derivative of a Stokeslet; all other
% tokens and their duplicate rules follow image_family_list.
if nargin==0
    self_test();
    return
end
if isempty(value)
    families = strings(1,0);
    return
end
value = string(value);
if any(ismissing(value),'all')
    error('ring_family_list:BadFamilies','Family tokens must not be missing.');
end
spec = regexprep(char(join(value(:).',"")),'\s','');
tokens = regexp(spec,'S_y|T_n2|T_[nt]|Q_s|[SRDT]','match');
if isempty(tokens) || ~strcmp(strjoin(tokens,''),spec)
    error('ring_family_list:BadFamilies','Invalid ring source family combination.');
end
families = string(tokens);
if nnz(families=="S_y")>1
    error('ring_family_list:BadFamilies','S_y may appear at most once.');
end
other = families(families~="S_y");
if ~isempty(other), image_family_list(other); end
end

function self_test()
%SELF_TEST Verify mixed token ordering and reject duplicates.
assert(isequal(ring_family_list("ST_nS_y"),["S","T_n","S_y"]));
assert(isequal(ring_family_list("S_y"),"S_y"));
assert(isempty(ring_family_list(strings(1,0))));
for invalid = ["S_yS_y","S_ySS","S_x"]
    try
        ring_family_list(invalid);
        error('ring_family_list:SelfTest','Invalid token accepted.');
    catch err
        assert(~strcmp(err.identifier,'ring_family_list:SelfTest'));
    end
end
fprintf('ring_family_list self-test passed.\n');
end
