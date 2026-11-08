function [figureHandle,report] = bubble_plot_collocation_residuals( ...
        solution,opt,visibility)
%BUBBLE_PLOT_COLLOCATION_RESIDUALS Show physical residuals at all fitted nodes.
% Include global nodes, graded rings, and source-depth patches. The four
% conditions use the validation plot's physical scales, before LS weights.
if nargin<3, visibility = 'on'; end
collocation = solution.grid;
collocation.count = collocation.collocationCount;
report = bubble_validate(solution,collocation,opt);
figureHandle = bubble_plot_residuals(report,collocation,solution,opt, ...
    visibility,'Collocation points');
end
