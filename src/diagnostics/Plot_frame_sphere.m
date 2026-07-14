function fig = Plot_frame_sphere(psi_a, psi_b, L_band, method, varargin)
%PLOT_FRAME_SPHERE
% Plot density and phase of two wavefunctions on the sphere.
%
% INPUTS
%   psi_a   : wavefunction of component a
%   psi_b   : wavefunction of component b
%   L_band  : spherical harmonic band-limit
%   method  : SSHT sampling method (e.g. 'MW')
%
% OPTIONAL NAME-VALUE PAIRS
%   'Visible'    : 'on' (default) or 'off'
%   'FigureTitle': title for the whole figure (default: '')
%
% OUTPUT
%   fig     : figure handle
%

p = inputParser;
addParameter(p, 'Visible', 'on', @(x) ischar(x) || isstring(x));
addParameter(p, 'FigureTitle', '', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});

fig_visible = char(p.Results.Visible);
fig_title   = char(p.Results.FigureTitle);

old_vis = get(groot, 'DefaultFigureVisible');
cleanupObj = onCleanup(@() set(groot, 'DefaultFigureVisible', old_vis));

if strcmpi(fig_visible, 'off')
    set(groot, 'DefaultFigureVisible', 'off');
end

fig = figure('Visible', fig_visible, 'Position', [100, 100, 1400, 900]);

tl = tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');

% -----------------------------
% Top-left: density of psi_a
% -----------------------------
ax1 = nexttile(tl, 1);
axes(ax1);
set(fig,'CurrentAxes',ax1)
ssht_plot_sphere(abs(psi_a).^2, L_band, ...
    'Method', method, ...
    'Type', 'colour', ...
    'ColourBar', false);
title('$|\psi_a|^2$', 'Interpreter', 'latex', 'FontSize', 18);
colormap(ax1, gray);
set(ax1, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
xlabel('$x/R$','Interpreter','latex','FontSize',15)
ylabel('$y/R$','Interpreter','latex','FontSize',15)
zlabel('$z/R$','Interpreter','latex','FontSize',15)

% -----------------------------
% Top-right: phase of psi_a
% -----------------------------
ax2 = nexttile(tl, 2);
axes(ax2);
set(fig,'CurrentAxes',ax2)
ssht_plot_sphere(angle(psi_a), L_band, ...
    'Method', method, ...
    'Type', 'colour', ...
    'ColourBar', false);
title('$\arg(\psi_a)$', 'Interpreter', 'latex', 'FontSize', 18);
colormap(ax2, hsv);
set(ax2, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
clim(ax2, [-pi pi]);
xlabel('$x/R$','Interpreter','latex','FontSize',15)
ylabel('$y/R$','Interpreter','latex','FontSize',15)
zlabel('$z/R$','Interpreter','latex','FontSize',15)

% -----------------------------
% Bottom-left: density of psi_b
% -----------------------------
ax3 = nexttile(tl, 3);
axes(ax3);
set(fig,'CurrentAxes',ax3)
ssht_plot_sphere(abs(psi_b).^2, L_band, ...
    'Method', method, ...
    'Type', 'colour', ...
    'ColourBar', false);
title('$|\psi_b|^2$', 'Interpreter', 'latex', 'FontSize', 18);
colormap(ax3, gray);
set(ax3, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
xlabel('$x/R$','Interpreter','latex','FontSize',15)
ylabel('$y/R$','Interpreter','latex','FontSize',15)
zlabel('$z/R$','Interpreter','latex','FontSize',15)

% -----------------------------
% Bottom-right: phase of psi_b
% -----------------------------
ax4 = nexttile(tl, 4);
axes(ax4);
set(fig,'CurrentAxes',ax4)
ssht_plot_sphere(angle(psi_b), L_band, ...
    'Method', method, ...
    'Type', 'colour', ...
    'ColourBar', false);
title('$\arg(\psi_b)$', 'Interpreter', 'latex', 'FontSize', 18);
colormap(ax4, hsv);
set(ax4, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
clim(ax4, [-pi pi]);
xlabel('$x/R$','Interpreter','latex','FontSize',15)
ylabel('$y/R$','Interpreter','latex','FontSize',15)
zlabel('$z/R$','Interpreter','latex','FontSize',15)

if ~isempty(fig_title)
    sgtitle(fig_title, 'Interpreter', 'latex', 'FontSize', 18);
end

set(fig, 'Visible', fig_visible);
drawnow
end