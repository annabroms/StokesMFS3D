function fig = ring_plot_residuals(points,residuals,labels,visibility)
%RING_PLOT_RESIDUALS Surface maps of normalized residuals at validation nodes.
%
%   fig = ring_plot_residuals(points,residuals,labels,visibility)
%
% Inputs:  points     Mv-by-3 Cartesian validation nodes, or a 1-by-K cell of
%                      Mk-by-3 node sets (one per residual, e.g. one region each);
%          residuals  1-by-K cell, each Mv-by-1 nonnegative, already normalized;
%          labels     1-by-K cell of titles stating each normalization;
%          visibility 'on' (default) or 'off'.
% Output:  fig, one row of K scatter3 maps colored by log10(max(residual,1e-16)).
if nargin == 0
    self_test();
    return
end
if nargin<4, visibility = 'on'; end
K = numel(residuals);
if ~iscell(points), points = repmat({points},1,numel(residuals)); end
fig = figure('Color','w','Visible',visibility,'Position',[40 80 420*K 420]);
layout = tiledlayout(fig,1,K,'TileSpacing','compact','Padding','compact');
for k = 1:K
    ax = nexttile(layout);
    scatter3(ax,points{k}(:,1),points{k}(:,2),points{k}(:,3),12, ...
        log10(max(residuals{k},1e-16)),'filled');
    axis(ax,'equal'); grid(ax,'on'); view(ax,35,22); colormap(ax,turbo);
    cb = colorbar(ax); cb.Label.String = 'log_{10} normalized residual';
    title(ax,labels{k},'Interpreter','none');
    xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
end
end

function self_test()
%SELF_TEST Build a hidden figure from random data.
p = randn(50,3); p = p./vecnorm(p,2,2);
fig = ring_plot_residuals(p,{abs(p(:,1)),abs(p(:,2))},{'a','b'},'off');
assert(numel(findobj(fig,'Type','Scatter'))==2);
close(fig);
fprintf('ring_plot_residuals self-test passed.\n');
end
