function families = image_family_list(value)
%IMAGE_FAMILY_LIST Parse ordered S/R/D/T/T_n/T_t source-family tokens.
if nargin == 0
    self_test(); families = strings(1,0); return
end
value = string(value);
if isempty(value) || any(ismissing(value))
    error('image_family_list:BadFamilies','Supply at least one source family.');
end
spec = regexprep(char(join(value(:).',"")),'\s','');
tokens = regexp(spec,'T_[nt]|[SRDT]','match');
families = string(tokens);
if isempty(tokens) || ~strcmp(strjoin(tokens,''),spec) || ...
        numel(unique(families))~=numel(families)
    error('image_family_list:BadFamilies', ...
        'Use distinct S, R, D, T, T_n, T_t tokens, e.g. "SDRT_nT_t".');
end
if any(families=="T") && any(ismember(families,["T_n","T_t"]))
    error('image_family_list:DuplicateStresslet', ...
        'Select either T or explicit T_n/T_t blocks to avoid duplicate directions.');
end
end

function self_test()
% Preserve compact legacy strings and distinguish the two stresslet tokens.
if ~isequal(image_family_list("SRDT"),["S","R","D","T"]) || ...
        ~isequal(image_family_list("SDRT_n T_t"),["S","D","R","T_n","T_t"]) || ...
        ~isequal(image_family_list(["S","T_n","T_t"]),["S","T_n","T_t"])
    error('image_family_list:SelfTest','Incorrect ordered family tokens.');
end
bad = ["SRSS","T_x","SD?","TT_n",""];
for value = bad
    rejected = false;
    try
        image_family_list(value);
    catch err
        rejected = startsWith(err.identifier,'image_family_list:');
    end
    if ~rejected, error('image_family_list:SelfTest','Invalid input accepted.'); end
end
fprintf('image_family_list self-test passed.\n');
end
