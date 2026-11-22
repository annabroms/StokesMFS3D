function export_readme_figures(outDir,selection)
%EXPORT_README_FIGURES Regenerate the illustrations shown in README.md.
%   EXPORT_README_FIGURES runs the demos unchanged, restyles their figures
%   with LaTeX-rendered ticks, labels, titles, colorbars and legends, and
%   exports PNG files to <repo>/figures.
%
%   EXPORT_README_FIGURES(outDir) writes to outDir instead.
%   EXPORT_README_FIGURES(outDir,selection) exports a subset; selection is a string
%   array with entries "ellipsoids", "wall" and/or "images".
%
%   Figures:
%     readme_ellipsoids.png     - ellipsoid_mobility_run(100,0.5,20,1,0,[],''),
%                                 redrawn at a finer plotting resolution
%     readme_wall_spheres.png   - demo_random_wall_torque (P = 40)
%     readme_image_cluster.png  - demo_cluster_resistance_with_images (P = 5)
%
%   The image demo is a script that calls "clear; close all", so it runs
%   first and in its own workspace.
%
% Anna Broms, Nov 22 2026

repoRoot = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || isempty(outDir)
    outDir = fullfile(repoRoot,'figures');
end
if nargin < 2
    selection = ["images","ellipsoids","wall"];
end
selection = string(selection);
if ~isfolder(outDir)
    mkdir(outDir);
end
run(fullfile(repoRoot,'setpath.m'));

if any(selection == "images")
    [fig,info] = run_image_cluster_demo();
    ax = findobj(fig,'Type','axes');
    set_latex(ax.Title,sprintf(['$P=%d$ spheres, gap $%s$, image-enhanced: ', ...
        'max BC error $%s$'],info.P,latex_number(info.gap), ...
        latex_number(info.maxBC)));
    set_latex(ax.Colorbar.Label,'$\log_{10}$ (BC error / max velocity)');
    label_xyz(ax);
    set(fig,'Position',[100,100,900,700]);
    apply_latex_style(fig,14);
    export(fig,fullfile(outDir,'readme_image_cluster.png'));
end

if any(selection == "ellipsoids")
    ellipsoid_mobility_run(100,0.5,20,1,0,[],'');
    oldFig = gcf; % the visualisation is the demo's last figure
    oldAx = findobj(oldFig,'Type','axes');
    fig = figure('Color','w');
    ax = axes(fig);
    hold(ax,'on');
    % Redraw each ellipsoid at a finer plotting resolution (visual only),
    % with lower ambient and stronger diffuse/specular light for depth.
    vertices = zeros(0,3);
    for s = findobj(oldAx,'Type','surface').'
        [X,Y,Z,C] = refine_ellipsoid_surface(s,48);
        surf(ax,X,Y,Z,C,'EdgeColor','none','FaceColor','interp', ...
            'FaceLighting','gouraud','AmbientStrength',0.35, ...
            'DiffuseStrength',0.7,'SpecularStrength',0.15, ...
            'SpecularExponent',25,'SpecularColorReflectance',0.5);
        vertices = [vertices;X(:),Y(:),Z(:)]; %#ok<AGROW>
    end
    close(oldFig);
    colormap(ax,parula);
    axis(ax,'equal','off');
    % Fit the camera to the cluster's projected extent in the default 3D
    % view and size the axes to that aspect ratio, so it is centred.
    [az,el] = deal(-37.5,30);
    T = viewmtx(az,el);
    screen = vertices*T(1:2,1:3).';
    lo = min(screen); hi = max(screen);
    halfExtent = (hi-lo)/2;
    centre = mean([min(vertices);max(vertices)]);
    radius = max(vecnorm(vertices-centre,2,2));
    depth = T(3,1:3);
    target = ((lo+hi)/2)*T(1:2,1:3) + (centre*depth.')*depth;
    distance = 10*radius;
    set(ax,'Projection','orthographic', ...
        'XLim',centre(1)+radius*[-1,1],'YLim',centre(2)+radius*[-1,1], ...
        'ZLim',centre(3)+radius*[-1,1],'CameraTarget',target, ...
        'CameraPosition',target+distance*depth,'CameraUpVector',[0,0,1], ...
        'CameraViewAngle',2*atand(1.03*min(halfExtent)/distance));
    height = 700;
    width = round(height*halfExtent(1)/halfExtent(2));
    set(fig,'Units','pixels','Position',[100,100,width+190,height+20]);
    set(ax,'Units','pixels','Position',[5,10,width,height]);
    % Key light from the upper left, dim fill light from the camera.
    keyLight = lightangle(ax,-60,45);
    keyLight.Color = 0.8*[1,1,1];
    light(ax,'Style','infinite','Position',depth,'Color',0.25*[1,1,1]);
    cb = colorbar(ax,'Units','pixels', ...
        'Position',[width+40,0.12*height,22,0.76*height]);
    set_latex(cb.Label,'Surface speed $|\mathbf{u}|$ from computed rigid-body motion');
    apply_latex_style(fig,14);
    export(fig,fullfile(outDir,'readme_ellipsoids.png'));
end

if any(selection == "wall")
    results = demo_random_wall_torque();
    close(results.geometryFigure);
    % Rebuild the surface-speed and residual panels in a two-panel figure;
    % the demo's three-panel layout is too narrow for the README.
    oldAxes = findobj(results.velocityFigure,'Type','axes');
    oldTitles = lower(string(arrayfun(@(a) string(a.Title.String),oldAxes)));
    speedAx = oldAxes(contains(oldTitles,'surface speed'));
    residualAx = oldAxes(contains(oldTitles,'residual'));
    fig = figure('Position',[100,100,1800,620],'Color','w');
    layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    panels = {speedAx,'Rigid-body surface speed','$|\mathbf{u}|$'; ...
        residualAx,sprintf('Independent surface residual, max $%s$', ...
            latex_number(results.details.validation.error)), ...
        '$\log_{10}$ relative residual'};
    for k = 1:2
        src = panels{k,1};
        ax = nexttile(layout);
        hold(ax,'on');
        copyobj(findobj(src.Children,'flat','-not','Type','light'),ax);
        colormap(ax,colormap(src));
        clim(ax,clim(src));
        axis(ax,'equal','tight');
        set(ax,'Projection','orthographic','Box','off', ...
            'XTick',-40:10:40,'YTick',-40:10:40,'ZTick',[0,5], ...
            'GridAlpha',0.12);
        grid(ax,'on');
        view(ax,-30,28);
        camlight(ax,'headlight');
        lighting(ax,'gouraud');
        cb = colorbar(ax);
        set_latex(cb.Label,panels{k,3});
        set_latex(ax.Title,panels{k,2});
        label_xyz(ax);
    end
    set_latex(layout.Title,sprintf(['$P=%d$ spheres above a no-slip wall, ', ...
        'each driven by a unit torque about $z$ (%d GMRES iterations)'], ...
        results.P,results.iterations));
    layout.Title.FontSize = 16;
    close(results.velocityFigure);
    apply_latex_style(fig,14);
    export(fig,fullfile(outDir,'readme_wall_spheres.png'));
end
end

% Resample one demo ellipsoid surface on a finer sphere(n) grid. The demo
% draws c + M*p for p on a coarse sphere grid, coloured by the rigid-body
% speed |v + w x (x - c)|, whose square is quadratic in p. Both are
% recovered exactly by least squares, so only the plotting changes.
function [X,Y,Z,C] = refine_ellipsoid_surface(s,n)
coarse = size(s.XData,1)-1;
[xs,ys,zs] = sphere(coarse);
p = [xs(:),ys(:),zs(:)];
affine = [p,ones(size(p,1),1)] \ [s.XData(:),s.YData(:),s.ZData(:)];
quad = @(p) [ones(size(p,1),1),p,p.^2,p(:,1).*p(:,2),p(:,1).*p(:,3),p(:,2).*p(:,3)];
speedCoeff = pinv(quad(p))*(double(s.CData(:)).^2);
[xf,yf,zf] = sphere(n);
pf = [xf(:),yf(:),zf(:)];
xyz = [pf,ones(size(pf,1),1)]*affine;
X = reshape(xyz(:,1),size(xf));
Y = reshape(xyz(:,2),size(xf));
Z = reshape(xyz(:,3),size(xf));
C = reshape(sqrt(max(quad(pf)*speedCoeff,0)),size(xf));
end

% Set a text object's interpreter before its string, to avoid TeX warnings.
function set_latex(textObject,str)
set(textObject,'Interpreter','latex','String',str);
end

function label_xyz(ax)
set_latex(ax.XLabel,'$x$');
set_latex(ax.YLabel,'$y$');
set_latex(ax.ZLabel,'$z$');
end

% Run the image demo script in this function's own workspace.
function [fig,info] = run_image_cluster_demo()
run('demo_cluster_resistance_with_images.m');
info = struct('P',P,'gap',gap,'maxBC',maxBC,'capMode',capMode); %#ok<NODEF>
end

% Render all figure text with the LaTeX interpreter, as in TikZ/pgfplots.
function apply_latex_style(fig,fontSize)
set(fig,'Color','w');
for ax = findall(fig,'Type','axes').'
    set(ax,'TickLabelInterpreter','latex','FontSize',fontSize);
    if ~isempty(ax.Toolbar)
        ax.Toolbar.Visible = 'off';
    end
    labels = [ax.Title,ax.Subtitle,ax.XLabel,ax.YLabel,ax.ZLabel];
    set(labels,'Interpreter','latex');
end
for cb = findall(fig,'Type','colorbar').'
    set(cb,'TickLabelInterpreter','latex','FontSize',fontSize);
    set(cb.Label,'Interpreter','latex','FontSize',fontSize);
end
for lg = findall(fig,'Type','legend').'
    set(lg,'Interpreter','latex','FontSize',fontSize);
end
for tl = findall(fig,'Type','tiledlayout').'
    set([tl.Title,tl.Subtitle,tl.XLabel,tl.YLabel],'Interpreter','latex');
end
drawnow;
end

function export(fig,file)
exportgraphics(fig,file,'Resolution',150,'BackgroundColor','white');
fprintf('Wrote %s\n',file);
end

% Format a number as LaTeX, e.g. 1.69\cdot 10^{-3}.
function s = latex_number(x)
e = floor(log10(abs(x)));
m = x/10^e;
if abs(m-1) < 5e-3
    s = sprintf('10^{%d}',e);
else
    s = sprintf('%.2f\\cdot 10^{%d}',m,e);
end
end
