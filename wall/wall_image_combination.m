function selected = wall_image_combination(combination)
%WALL_IMAGE_COMBINATION Decode a readable image-source combination string.
% Order: [S,R,T_n,T_n2,T_t,D]. T expands to all three stresslet blocks.
% Examples: "SRTD", "ST_tD", "ST_nT_n2D", "T_nT_n2", "R+T_t".
% T_n and T_n2 are normal to the line (tangent to the wall); T_t is along
% the line (normal to the wall). Standard proxy S is implicit in every solve.
% No-argument invocation tests full, split and invalid specifications.

if nargin == 0
    self_test(); selected = []; return
end
combination = string(combination);
if ~isscalar(combination) || ismissing(combination) || strlength(combination) == 0
    error('wall:ImageCombination','Supply one nonempty combination string.');
end
remaining = erase(combination,["+"," ",","]);
selected = false(1,6);
% Consume longer tokens before their prefixes: T_n2 must precede T_n.
tokens = ["T_n2","T_n","T_t","T","S","R","D"];
indices = [4,3,5,0,1,2,6];
for k = 1:numel(tokens)
    if contains(remaining,tokens(k))
        if indices(k) == 0
            selected(3:5) = true;
        else
            selected(indices(k)) = true;
        end
        remaining = erase(remaining,tokens(k));
    end
end
if strlength(remaining) > 0 || ~any(selected)
    error('wall:ImageCombination','Unknown source token in "%s".',combination);
end
end

function self_test()
% Split orientations must not accidentally select every token containing T.
names = ["SRTD","T_t","T_nT_n2","ST_tD","R+T_n+T_n2"];
expected = logical([1,1,1,1,1,1;0,0,0,0,1,0;0,0,1,1,0,0; ...
    1,0,0,0,1,1;0,1,1,1,0,0]);
for k = 1:numel(names)
    if ~isequal(wall_image_combination(names(k)),expected(k,:))
        error('wall:ImageCombinationTest','Combination parsing failed for %s.',names(k));
    end
end
caught = false;
try
    wall_image_combination("ST_bad");
catch exception
    caught = strcmp(exception.identifier,'wall:ImageCombination');
end
if ~caught, error('wall:ImageCombinationTest','Invalid tokens were accepted.'); end
fprintf('wall_image_combination PASS: full and split stresslet selections.\n');
end
