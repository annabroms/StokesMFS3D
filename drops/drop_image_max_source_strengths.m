function [innerProxy,outerProxy,innerImages,outerImages] = ...
    drop_image_max_source_strengths(grid,lambda,gamma,families)
%DROP_IMAGE_MAX_SOURCE_STRENGTHS Max physical vector norm by source family.
if nargin == 0
    self_test();
    innerProxy = []; outerProxy = []; innerImages = []; outerImages = [];
    return
end
families = image_family_list(families);
innerProxy = 0;
outerProxy = 0;
innerImages = nan(1,numel(families));
outerImages = nan(1,numel(families));

% The exterior-flow basis is inside each drop; the interior-flow basis is outside.
for i = 1:numel(grid.body)
    bases = {grid.body(i).basis,grid.body(i).interiorBasis};
    strengths = {lambda{i},gamma{i}};
    for side = 1:2
        blocks = bases{side};
        coeff = image_source_coefficients(blocks,strengths{side});
        for b = 1:numel(blocks)
            if isempty(coeff{b}), continue; end
            value = max(vecnorm(coeff{b},2,2));
            if b==1
                if side==1
                    innerProxy = max(innerProxy,value);
                else
                    outerProxy = max(outerProxy,value);
                end
                continue
            end
            family = find(families==blocks(b).family,1);
            if isempty(family)
                error('drop_image_max_source_strengths:UnknownFamily', ...
                    'A basis family is absent from the requested list.');
            end
            if side==1
                if isnan(innerImages(family)) || value>innerImages(family)
                    innerImages(family) = value;
                end
            else
                if isnan(outerImages(family)) || value>outerImages(family)
                    outerImages(family) = value;
                end
            end
        end
    end
end
end

function self_test()
% Check interleaving, basis scaling, and absent outer image families.
b = struct('family',"S",'points',[0 0 0], ...
    'normalDirections',zeros(0,3),'scale',1);
body.basis = repmat(b,3,1);
body.interiorBasis = repmat(b,3,1);
body.basis(2).family = "S";
body.basis(3).family = "R";
body.interiorBasis(2).family = "S";
body.interiorBasis(3).family = "R";
body.interiorBasis(3).scale = 2;
grid.body = body;
lambda = {[3;4;0;0;0;2;1;1;1]};
gamma = {[0;0;6;0;0;7;1;0;0]};
[ip,op,ii,oi] = drop_image_max_source_strengths( ...
    grid,lambda,gamma,"SR");
if abs(ip-5)>1e-14 || abs(op-6)>1e-14 || ...
        norm(ii-[2 sqrt(3)])>1e-14 || norm(oi-[7 2])>1e-14
    error('drop_image_max_source_strengths:SelfTest', ...
        'The source-family maxima do not match the known strengths.');
end
grid.body.interiorBasis = b;
gamma = {[0;0;6]};
[~,~,~,oi] = drop_image_max_source_strengths(grid,lambda,gamma,"SR");
if ~all(isnan(oi))
    error('drop_image_max_source_strengths:SelfTest', ...
        'Absent outer image families must be reported as NaN.');
end
grid.body.basis(2).family = "T_n";
grid.body.basis(3).family = "T_t";
[~,~,ii] = drop_image_max_source_strengths(grid,lambda,gamma,"SDRT_nT_t");
if ~all(isnan(ii(1:3))) || norm(ii(4:5)-[2 sqrt(3)])>1e-14
    error('drop_image_max_source_strengths:SelfTest','Stresslet tokens were not separated.');
end
fprintf('drop_image_max_source_strengths self-test passed.\n');
end
