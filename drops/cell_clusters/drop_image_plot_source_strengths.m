function fig = drop_image_plot_source_strengths(grid,lambda,gamma,label,visibility)
%DROP_IMAGE_PLOT_SOURCE_STRENGTHS Plot log10 strength magnitude at each source.
if nargin == 0
    self_test(); fig = []; return
end
if nargin < 4, label = ''; end
if nargin < 5, visibility = 'on'; end

% Each source block has three interleaved coefficients per point.
families = strings(1,0);
for i = 1:numel(grid.body)
    for b = 2:numel(grid.body(i).basis)
        families(end+1) = grid.body(i).basis(b).family; %#ok<AGROW>
    end
    for b = 2:numel(grid.body(i).interiorBasis)
        families(end+1) = grid.body(i).interiorBasis(b).family; %#ok<AGROW>
    end
end
families = unique(families,'stable');
hasOuterImages = any(arrayfun( ...
    @(b) numel(b.interiorBasis)>1,grid.body));
fig = figure('Name','Drop source strengths','Color','w', ...
    'Visible',visibility);
if hasOuterImages
    nCols = 1+numel(families);
    nPanels = 2*nCols;
    tiledlayout(fig,2,nCols,'TileSpacing','compact');
else
    nPanels = 2+numel(families);
    nCols = ceil(sqrt(nPanels));
    tiledlayout(fig,ceil(nPanels/nCols),nCols,'TileSpacing','compact');
end

% Separate source ownership so mirrored strengths remain visible.
for panel = 1:nPanels
    if hasOuterImages
        side = 1+floor((panel-1)/nCols);
        familyIndex = mod(panel-1,nCols);
    elseif panel==1
        side = 1; familyIndex = 0;
    elseif panel==2
        side = 2; familyIndex = 0;
    else
        side = 1; familyIndex = panel-2;
    end
    ax = nexttile;
    hold(ax,'on');
    for i = 1:numel(grid.body)
        if side==1
            basis = grid.body(i).basis;
            strength = lambda{i};
            marker = 'o';
        else
            basis = grid.body(i).interiorBasis;
            strength = gamma{i};
            marker = '^';
        end
        coeff = image_source_coefficients(basis,strength);
        for b = 1:numel(basis)
            if familyIndex==0 && b~=1 || familyIndex>0 && ...
                    (b==1 || basis(b).family~=families(familyIndex))
                continue
            end
            x = basis(b).points;
            if isempty(x), continue; end
            magnitude = vecnorm(coeff{b},2,2);
            value = log10(max(magnitude,realmin));
            scatter3(ax,x(:,1),x(:,2),x(:,3), ...
                20+20*(familyIndex>0),value,'filled','Marker',marker);
        end
    end
    axis(ax,'equal'); view(ax,3);
    set(ax,'XGrid','on','YGrid','on','ZGrid','on');
    xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
    if familyIndex==0 && side==1
        title(ax,'Exterior-flow proxy');
    elseif familyIndex==0
        title(ax,'Interior-flow proxy');
    elseif side==1
        title(ax,sprintf('%s interior images',families(familyIndex)));
    else
        title(ax,sprintf('%s exterior images',families(familyIndex)));
    end
    cb = colorbar(ax);
    cb.Label.String = 'log_{10} ||strength||_2';
end
colormap(fig,turbo);
sgtitle(fig,sprintf('%s source strengths',label));
drawnow;
end

function self_test()
% Check every plotted node and a known three-component norm.
q = [0 0 0;2.1 0 0]; a = ones(2,1);
opt = drop_image_options(struct('N',30,'N_image',1, ...
    'imageFamilies',"SR",'capPointCounts',[12 12], ...
    'mirror_ext',true),2,a);
grid = drop_image_discretization(q,a,opt);
lambda = cell(2,1); gamma = cell(2,1);
expectedCount = 0;
for i = 1:2
    ni = sum(arrayfun(@(b) size(b.points,1),grid.body(i).basis));
    no = sum(arrayfun(@(b) size(b.points,1),grid.body(i).interiorBasis));
    lambda{i} = zeros(3*ni,1);
    gamma{i} = zeros(3*no,1);
    expectedCount = expectedCount+ni+no;
end
lambda{1}(1:3) = [3;4;0];
outerStart = 3*size(grid.body(1).interiorBasis(1).points,1);
gamma{1}(outerStart+(1:3)) = [0;0;7];
fig = drop_image_plot_source_strengths( ...
    grid,lambda,gamma,'self-test','off');
plots = findall(fig,'Type','Scatter');
actualCount = sum(arrayfun(@(h) numel(h.XData),plots));
values = arrayfun(@(h) h.CData(:),plots,'UniformOutput',false);
allValues = vertcat(values{:});
if actualCount~=expectedCount || ...
        numel(findall(fig,'Type','Axes'))~=6 || ...
        ~any(abs(allValues-log10(5))<1e-12) || ...
        ~any(abs(allValues-log10(7))<1e-12)
    error('drop_image_plot_source_strengths:SelfTest', ...
        'Source positions or vector magnitudes were plotted incorrectly.');
end
axesList = findall(fig,'Type','Axes');
titles = arrayfun(@(ax) string(ax.Title.String),axesList);
outerAxis = axesList(titles=="S exterior images");
outerPlots = findall(outerAxis,'Type','Scatter');
outerValues = arrayfun(@(h) h.CData(:),outerPlots,'UniformOutput',false);
if numel(outerAxis)~=1 || ...
        ~any(abs(vertcat(outerValues{:})-log10(7))<1e-12)
    error('drop_image_plot_source_strengths:SelfTest', ...
        'Mirrored image strengths are missing from their panel.');
end
close(fig);
% Each explicit stresslet direction needs its own panel on both sides.
for i = 1:2
    grid.body(i).basis(2).family = "T_n";
    grid.body(i).basis(3).family = "T_t";
    grid.body(i).interiorBasis(2).family = "T_n";
    grid.body(i).interiorBasis(3).family = "T_t";
end
fig = drop_image_plot_source_strengths(grid,lambda,gamma,'joint','off');
axesList = findall(fig,'Type','Axes');
titles = arrayfun(@(ax) string(ax.Title.String),axesList);
expectedTitles = ["T_n interior images","T_t interior images", ...
    "T_n exterior images","T_t exterior images"];
if ~all(ismember(expectedTitles,titles))
    error('drop_image_plot_source_strengths:SelfTest','Missing stresslet direction panel.');
end
close(fig);
opt.useImages = false;
grid = drop_image_discretization(q,a,opt);
for i = 1:2
    lambda{i} = zeros(3*size(grid.body(i).basis(1).points,1),1);
    gamma{i} = zeros(3*size(grid.body(i).interiorBasis(1).points,1),1);
end
fig = drop_image_plot_source_strengths( ...
    grid,lambda,gamma,'no images','off');
if numel(findall(fig,'Type','Axes'))~=2
    error('drop_image_plot_source_strengths:SelfTest', ...
        'The no-image case should show only the two proxy panels.');
end
close(fig);
fprintf('drop_image_plot_source_strengths self-test passed.\n');
end
