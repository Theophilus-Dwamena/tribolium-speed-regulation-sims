%% speed_regulation_three_gradients.m
%
% Simulates embryonic gene expression for a SINGLE GENE driven purely by
% a speed gradient.
%
% The morphogen gradient modulates how fast each cell accumulates phase y(x,t).
% When y crosses successive thresholds it changes gene state (0 to 1).
%
% Switch between gradient forms by setting:
%   gradientForm = 1   --> gradient-based      (static linear gradient, constant in time)
%   gradientForm = 2   --> graded wavefront     (smooth retracting step)
%   gradientForm = 3   --> wavefront-based      (sharp retracting step)
%
% Outputs:
%   - Animated movie:  speed_regulation_<gradientLabel>_animation.mp4
%   - Kymograph PDF:   speed_regulation_<gradientLabel>_kymographs.pdf
%   - Snapshot PDF:    speed_regulation_<gradientLabel>_snapshot_t<t_snapshot>.pdf

clear; clc; close all;

realtime_factor = 1;

%% ── Spatial and temporal domains ─────────────────────────────────────────────
x_min = 0;
x_max = 1;
Nx    = 50000;
x     = linspace(x_min, x_max, Nx);

t_min = 0;
t_max = 5;
dt    = 0.01;
time  = t_min:dt:t_max;
Nt    = length(time);

%% ── Gene state parameters ────────────────────────────────────────────────────
N_states = 1;
lambda   = 1;      % exon decay / processing rate
y_on     = 0.3;    % phase at which gene switches ON
y_off    = 0.5;    % phase at which gene switches OFF

%% ── GLOBAL COLOR SCHEME ──────────────────────────────────────────────────────
C = struct();

C.bg     = [1 1 1];
C.off    = [1 1 1];

C.intron = [0.1 0.1 0.85];   % blue
C.exon   = [0.1 0.1 0.85];   % blue

C.grad   = [0.6 0.6 0.6];    % gray

C.stateMap = [C.off; C.intron];
C.traceMap = [C.off; C.exon  ];

%% ── Gradient form ────────────────────────────────────────────────────────────
%   1 = gradient-based   (static linear gradient, constant in time)
%   2 = graded wavefront (smooth retracting step)
%   3 = wavefront-based  (sharp retracting step)
gradientForm = 3;

% Human-readable label used in all figure titles and file names
switch gradientForm
    case 1,    gradientLabel = 'gradient-based';
    case 2,    gradientLabel = 'graded-wavefront';
    case 3,    gradientLabel = 'wavefront-based';
    otherwise, gradientLabel = 'unknown';
end

%% ── Gradient amplitude parameters ───────────────────────────────────────────
param_max = 1;
param_min = 0.0;
gradStart = 0;
gradEnd   = x_max;
v         = 1.5 * 0.12;   % retraction speed [space / time]
stepWidth = 0.0025;        % spatial width of the linear transition zone

%% ── Speed (phase accumulation rate) ─────────────────────────────────────────
default_c0 = 0.3;

%% ── Initial condition ────────────────────────────────────────────────────────
y    = zeros(Nx, 1);   % all cells start with phase = 0
exon = zeros(Nx, 1);   % trailing exon signal

%% ── Snapshot time ────────────────────────────────────────────────────────────
t_snapshot  = 2;
snapshot_ti = round((t_snapshot - t_min) / dt) + 1;
snap_gradVals  = zeros(1, Nx);
snap_geneState = zeros(Nx, 1);
snap_exon      = zeros(Nx, 1);

%% ── Kymograph storage ────────────────────────────────────────────────────────
gradientMatrix  = zeros(Nt, Nx);
geneStateMatrix = zeros(Nt, Nx);
exonMatrix      = zeros(Nt, Nx);

%% ── Animated figure setup ────────────────────────────────────────────────────
fig_anim = figure('Color', 'w', ...
    'Name', sprintf('Speed regulation — %s — spatial profiles', gradientLabel));

tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile([2 1]);

hGradient = plot(x, zeros(size(x)), '-',  'Color', C.grad,   'LineWidth', 2);
hold on;
hState    = plot(x, zeros(1, Nx),   '--', 'Color', C.intron, 'LineWidth', 2);
hExon     = plot(x, zeros(1, Nx),   '-',  'Color', C.exon,   'LineWidth', 2);
hold off;

lg = legend({'gradient', 'intron', 'exon'}, 'Location', 'northwest');
legend('boxoff');
lg.FontSize = 14;

ylim([-0.05 1.05]);
xlim([x_min x_max]);
xlabel('Spatial position');
ylabel('Value');
title(ax1, sprintf('Speed regulation — %s — spatial profiles', gradientLabel));
grid on;

%% ── Movie writer setup ───────────────────────────────────────────────────────
frame_interval = 10;
nFrames        = floor(Nt / frame_interval);
frames(nFrames) = struct('cdata', [], 'colormap', []);
frameCounter   = 1;

movie_file = sprintf('speed_regulation_%s_animation.mp4', gradientLabel);
vw = VideoWriter(movie_file, 'MPEG-4');
vw.FrameRate = max(1, round(1 / (dt * frame_interval)));
open(vw);

%% ── Time integration (Euler) ─────────────────────────────────────────────────
for ti = 1:Nt
    t_current = time(ti);
    gradVals  = zeros(1, Nx);

    for j = 1:Nx
        param_val = getGradientValue( ...
            x(j), t_current, gradientForm, ...
            param_max, param_min, gradStart, gradEnd, ...
            v, stepWidth);

        gradVals(j) = param_val;

        c_val = param_val * default_c0;

        if y(j) == 0
            newPhase = y(j) + dt * c_val;
            if newPhase > 0
                y(j) = newPhase;
            end
        else
            y(j) = y(j) + dt * c_val;
        end
    end

    gradientMatrix(ti, :) = gradVals;

    gene_state = zeros(Nx, 1);
    for j = 1:Nx
        gene_state(j) = findGeneState(y(j), N_states, y_on, y_off);
        x_intronic = gene_state(j);
        exon(j) = exon(j) + dt * (x_intronic - lambda * exon(j));
    end
    geneStateMatrix(ti, :) = gene_state;
    exonMatrix(ti, :)      = exon;

    if ti == snapshot_ti
        snap_gradVals  = gradVals;
        snap_geneState = gene_state;
        snap_exon      = exon;
    end

    if mod(ti, frame_interval) == 0
        set(hGradient, 'YData', gradVals);
        set(hState,    'YData', gene_state');
        set(hExon,     'YData', exon');
        title(ax1, sprintf('Speed regulation — %s — spatial profiles,  t = %.2f', ...
            gradientLabel, t_current));
        drawnow;
        pause(dt * frame_interval * realtime_factor);
        frames(frameCounter) = getframe(fig_anim);
        frameCounter = frameCounter + 1;
    end
end

%% ── Save movie ───────────────────────────────────────────────────────────────
for k = 1:length(frames)
    if ~isempty(frames(k).cdata)
        writeVideo(vw, frames(k));
    end
end
close(vw);
fprintf('Movie saved: %s\n', movie_file);

%% ── Kymograph figure ─────────────────────────────────────────────────────────
fig_kymo = figure('Color', 'w', ...
    'Name', sprintf('Speed regulation — %s — kymographs', gradientLabel));

ax1k = subplot(3, 1, 1);
imagesc(x, time, gradientMatrix);
set(ax1k, 'YDir', 'normal');
colormap(ax1k, flipud(gray));
colorbar(ax1k);
xlabel(ax1k, 'AP position');
ylabel(ax1k, 'Time');
title(ax1k, sprintf('Gradient  (%s)', gradientLabel));

ax2k = subplot(3, 1, 2);
imagesc(x, time, geneStateMatrix);
set(ax2k, 'YDir', 'normal');
colormap(ax2k, C.stateMap);
caxis(ax2k, [0 1]);
colorbar(ax2k);
xlabel(ax2k, 'AP position');
ylabel(ax2k, 'Time');
title(ax2k, sprintf('Intron  (%s)', gradientLabel));

ax3k = subplot(3, 1, 3);
imagesc(x, time, exonMatrix);
set(ax3k, 'YDir', 'normal');
colormap(ax3k, C.traceMap);
colorbar(ax3k);
xlabel(ax3k, 'AP position');
ylabel(ax3k, 'Time');
title(ax3k, sprintf('Exon  (%s)', gradientLabel));

kymo_file = sprintf('speed_regulation_%s_kymographs.pdf', gradientLabel);
exportgraphics(fig_kymo, kymo_file, 'ContentType', 'vector');
fprintf('Figure saved: %s\n', kymo_file);

%% ── Snapshot figure ──────────────────────────────────────────────────────────
fig_snap = figure('Color', 'w', ...
    'Name', sprintf('Speed regulation — %s — snapshot t = %d', gradientLabel, t_snapshot));

plot(x, snap_gradVals,  '-',  'Color', C.grad,   'LineWidth', 2); hold on;
plot(x, snap_geneState, '--', 'Color', C.intron,  'LineWidth', 2);
plot(x, snap_exon,      '-',  'Color', C.exon,    'LineWidth', 2);
hold off;

lg = legend({'gradient', 'intron', 'exon'}, 'Location', 'northwest');
legend('boxoff');
lg.FontSize = 14;
ylim([-0.05 1.05]);
xlim([x_min x_max]);
xlabel('Spatial position');
ylabel('Value');
title(sprintf('Speed regulation — %s — snapshot,  t = %d', gradientLabel, t_snapshot));
grid on;

snap_file = sprintf('speed_regulation_%s_snapshot_t%d.pdf', gradientLabel, t_snapshot);
exportgraphics(fig_snap, snap_file, 'ContentType', 'vector');
fprintf('Figure saved: %s\n', snap_file);

%% ════════════════════════════════════════════════════════════════════════════
%%  HELPER FUNCTIONS
%% ════════════════════════════════════════════════════════════════════════════

function val = getGradientValue(x, t, gradientForm, ...
        param_max, param_min, gradStart, gradEnd, v, stepWidth)

    switch gradientForm

        case 1  % gradient-based: static linear ramp
            if x <= gradStart
                val = param_min;
            elseif x >= gradEnd
                val = param_max;
            else
                val = param_min + (param_max - param_min) * ...
                      (x - gradStart) / (gradEnd - gradStart);
            end

        case 2  % graded wavefront: smooth retracting step
            step_center = gradStart + v * t;
            if stepWidth <= 0
                if x < step_center
                    val = param_min;
                else
                    val = param_max;
                end
            else
                left  = step_center - stepWidth / 2;
                right = step_center + stepWidth / 2;
                if x < left
                    val = param_min;
                elseif x > right
                    val = param_max;
                else
                    val = param_min + (param_max - param_min) * ...
                          (x - left) / stepWidth;
                end
            end

        case 3  % wavefront-based: sharp retracting step
            step_center = gradStart + v * t;
            if x < step_center
                val = param_min;
            else
                val = param_max;
            end

        otherwise
            val = param_max;
    end
end

function state = findGeneState(y_val, N, y_on, y_off)
    sigma  = (y_off - y_on) / 4;
    y_star = 0.5 * (y_on + y_off);
    state  = exp(-((y_val - y_star)^2) / (2 * sigma^2));
    state  = max(0, min(1, state));
end