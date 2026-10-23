function fig = plot_drop_image_source_family_comparison(results,outputFile)
%PLOT_DROP_IMAGE_SOURCE_FAMILY_COMPARISON Plot max physical source magnitudes.
% One marker per selected basis and source family; inner and mirrored outer
% lines have separate panels. Missing family tokens are left blank.
if nargin<2, outputFile = ""; end
families = ["S","D","R","T_n","T_t","T_n2","Q_s"];
labels = ["proxy S",families];
selected = [results.region]==0 & string({results.location})=="surface";
rows = results(selected);
requested = unique([rows.NRequested],'stable');
cases = unique(string({rows.family}),'stable');
colors = lines(numel(cases));
fig = figure('Name','Two-drop image source magnitudes','Color','w', ...
    'Position',[80 80 1800 1050]);
tiledlayout(fig,2,numel(requested),'TileSpacing','compact');
legendHandles = gobjects(numel(cases),1);
for side = 1:2
    for ni = 1:numel(requested)
        ax = nexttile((side-1)*numel(requested)+ni);
        hold(ax,'on');
        for k = 1:numel(cases)
            row = rows([rows.NRequested]==requested(ni) & ...
                string({rows.family})==cases(k));
            assert(isscalar(row));
            values = nan(1,numel(labels));
            if side==1
                values(1) = row.innerProxyMax;
                imageMax = row.innerImageMax;
            else
                values(1) = row.outerProxyMax;
                imageMax = row.outerImageMax;
            end
            tokens = image_family_list(row.family);
            for j = 1:numel(tokens)
                index = find(families==tokens(j),1);
                values(index+1) = imageMax(j);
            end
            x = (1:numel(labels))+0.11*(k-(numel(cases)+1)/2);
            h = semilogy(ax,x,values,'o','Color',colors(k,:), ...
                'MarkerFaceColor',colors(k,:),'MarkerSize',6, ...
                'DisplayName',char(cases(k)));
            if side==1 && ni==1, legendHandles(k) = h; end
        end
        xlim(ax,[0.5 numel(labels)+0.5]);
        set(ax,'YScale','log');
        xticks(ax,1:numel(labels)); xticklabels(ax,labels);
        xtickangle(ax,35); grid(ax,'on');
        ylabel(ax,'Maximum physical strength magnitude');
        if side==1, sourceSide = 'inner'; else, sourceSide = 'mirrored outer'; end
        title(ax,sprintf('%s, requested N=%d',sourceSide,requested(ni)));
        if side==1 && ni==1, firstAxes = ax; end
    end
end
key = legend(firstAxes,legendHandles,cellstr(cases), ...
    'Orientation','horizontal','Interpreter','none');
key.Layout.Tile = 'north';
if strlength(string(outputFile))>0
    outputFile = string(outputFile);
    folder = fileparts(outputFile);
    if strlength(folder)>0 && ~isfolder(folder), mkdir(folder); end
    exportgraphics(fig,outputFile,'Resolution',200);
end
end
