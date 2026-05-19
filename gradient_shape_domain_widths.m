function gradient_shape_domain_widths()
% =========================================================================
% gradient_shape_domain_widths.m
% -------------------------------------------------------------------------
% Single-gene speed-gradient simulation (no GRN, no exon).
% Dynamical system (from GeneSpeedGradient.m):
%   d(phi)/dt = G(x, t) * c0
%   expression = findGeneState(phi)   [Gaussian pulse in phase space]
%
% Two gradient shapes are compared simultaneously (gradientForm = 1 only):
%   Linear  gradient → solid  lines
%   Concave gradient → dashed lines
%
% Color scheme:
%   Gradient (linear)  → solid  gray   [0.60 0.60 0.60]
%   Gradient (concave) → dashed gray   [0.60 0.60 0.60]
%   Pulse    (linear)  → solid  black  [0    0    0   ]
%   Pulse    (concave) → dashed black  [0    0    0   ]
%
% Figures produced:
%   1. Animated spatial profiles  (2 gradient curves + 2 pulse curves)
%      saved as vector PDF: gradient_shape_domain_widths_spatial_profiles.pdf
%
% Movie saved as: gradient_shape_domain_widths_animation.mp4
% =========================================================================

clearvars; clc; close all;

realtime_factor = 1;    % set > 1 to slow animation down, < 1 to speed up

% ── 1) SPATIAL & TEMPORAL DOMAIN ─────────────────────────────────────────
x_min = 0;
x_max = 1;
Nx    = 5000;
x     = linspace(x_min, x_max, Nx);

t_min = 0;
t_max = 5;
dt    = 0.01;
time  = t_min : dt : t_max;
Nt    = length(time);

% ── DISPLAY SCALE ─────────────────────────────────────────────────────────
AP_display_scale = 1;
x_display  = x * AP_display_scale;
x_disp_min = x_min * AP_display_scale;
x_disp_max = x_max * AP_display_scale;

% ── 2) GRADIENT FORM ──────────────────────────────────────────────────────
gradientForm = 1;   % Static linear ramp vs concave power-law

% ── 3) GRADIENT SHAPE PARAMETERS ─────────────────────────────────────────
param_max      = 1.0;
param_min      = 0.0;
gradStart      = 0;
gradEnd        = 1;
v_wave         = 1.5 * 0.12;
retract_offset = 0.0;
stepWidth      = 0.5;
concave_p      = 0.5;
thresh_frac    = 0.5;
y_thresh       = 0.3;
b2g_time       = 0.0;

% ── 4) DYNAMICAL SYSTEM PARAMETERS ───────────────────────────────────────
c0    = 0.3;
y_on  = 0.3;
y_off = 0.5;

% ── 5) COLOR SCHEME ──────────────────────────────────────────────────────
gray_col  = [0.60 0.60 0.60];
black_col = [0    0    0   ];

% ── 6) GRADIENT LABELS ───────────────────────────────────────────────────
lin_label = 'static linear';
con_label = sprintf('static concave (p=%.2f)', concave_p);

% =========================================================================
% ── 7) EULER INTEGRATION ─────────────────────────────────────────────────
% =========================================================================
fprintf('Running Euler integration (%d steps × %d positions) ...\n', Nt, Nx);

gradMatrix_lin = zeros(Nt, Nx);
gradMatrix_con = zeros(Nt, Nx);
pulse_lin      = zeros(Nx, Nt);
pulse_con      = zeros(Nx, Nt);

phi_lin = zeros(Nx, 1);
phi_con = zeros(Nx, 1);

for ti = 1:Nt
    t_now = time(ti);

    G_lin_vec = arrayfun(@(xj) getGradientLinear(xj, t_now, gradientForm, ...
        param_max, param_min, gradStart, gradEnd, v_wave, retract_offset, ...
        stepWidth, b2g_time), x(:));

    G_con_vec = arrayfun(@(xj) getGradientConcave(xj, t_now, gradientForm, ...
        param_max, param_min, gradStart, gradEnd, v_wave, retract_offset, ...
        concave_p, stepWidth, thresh_frac, y_thresh, b2g_time), x(:));

    gradMatrix_lin(ti, :) = G_lin_vec';
    gradMatrix_con(ti, :) = G_con_vec';

    for j = 1:Nx
        dy_lin = G_lin_vec(j) * c0;
        if phi_lin(j) == 0
            np = phi_lin(j) + dt * dy_lin;
            if np > 0, phi_lin(j) = np; end
        else
            phi_lin(j) = phi_lin(j) + dt * dy_lin;
        end

        dy_con = G_con_vec(j) * c0;
        if phi_con(j) == 0
            np = phi_con(j) + dt * dy_con;
            if np > 0, phi_con(j) = np; end
        else
            phi_con(j) = phi_con(j) + dt * dy_con;
        end
    end

    pulse_lin(:, ti) = arrayfun(@(ph) findGeneState(ph, y_on, y_off), phi_lin);
    pulse_con(:, ti) = arrayfun(@(ph) findGeneState(ph, y_on, y_off), phi_con);
end

fprintf('Integration complete.\n');

% =========================================================================
% FIGURE 1 – Animated spatial profiles
% =========================================================================
fig1 = figure('Name', 'Gradient shape comparison — spatial profiles', 'Color', 'w');

hGrad_lin = plot(x_display, gradMatrix_lin(1,:), '-',  'Color', gray_col,  'LineWidth', 1.5);
hold on;
hGrad_con = plot(x_display, gradMatrix_con(1,:), '--', 'Color', gray_col,  'LineWidth', 1.5);
hPuls_lin = plot(x_display, pulse_lin(:,1)',     '-',  'Color', black_col, 'LineWidth', 2.5);
hPuls_con = plot(x_display, pulse_con(:,1)',     '--', 'Color', black_col, 'LineWidth', 2.5);
hold off;

legend({ ...
    ['gradient (' lin_label ')'], ...
    ['gradient (' con_label ')'], ...
    'pulse (linear)', ...
    'pulse (concave)'}, ...
    'Location', 'northwest');
legend('boxoff');

y_ceil = max([param_max, pulse_lin(:)', pulse_con(:)']) * 1.1;
ylim([-0.05, max(y_ceil, 0.1)]);
xlim([x_disp_min x_disp_max]);
xlabel('A-P position');
ylabel('Expression / Gradient value');
title('Gradient shape comparison — spatial profiles');
grid on; box off;

% ── Movie writer setup ────────────────────────────────────────────────────
frame_interval = 5;
nFrames        = floor(Nt / frame_interval);
frames(nFrames) = struct('cdata', [], 'colormap', []);
frameCounter   = 1;

outputFile = 'gradient_shape_domain_widths_animation.mp4';
movieFPS   = max(1, round(1 / (dt * frame_interval)));
vw = VideoWriter(outputFile, 'MPEG-4');
vw.FrameRate = movieFPS;
open(vw);

for ti = 1:Nt
    if mod(ti, frame_interval) == 0
        t_now = time(ti);
        set(hGrad_lin, 'YData', gradMatrix_lin(ti, :));
        set(hGrad_con, 'YData', gradMatrix_con(ti, :));
        set(hPuls_lin, 'YData', pulse_lin(:, ti)');
        set(hPuls_con, 'YData', pulse_con(:, ti)');
        title(sprintf('Gradient shape comparison — spatial profiles,  t = %.2f', t_now));
        drawnow;
        pause(dt * frame_interval * realtime_factor);
        if frameCounter <= nFrames
            frames(frameCounter) = getframe(fig1);
            writeVideo(vw, frames(frameCounter));
            frameCounter = frameCounter + 1;
        end
    end
end

close(vw);
fprintf('Movie saved: %s\n', outputFile);

% Save Figure 1 as vector PDF (final state of animation)
exportgraphics(fig1, 'gradient_shape_domain_widths_spatial_profiles.pdf', 'ContentType', 'vector');
fprintf('Figure saved: gradient_shape_domain_widths_spatial_profiles.pdf\n');

% =========================================================================
%                         HELPER FUNCTIONS
% =========================================================================

    function state = findGeneState(y_val, y_on, y_off)
        sigma  = (y_off - y_on) / 4;
        y_star = 0.5 * (y_on + y_off);
        state  = exp(-((y_val - y_star)^2) / (2 * sigma^2));
        state  = max(0, min(1, state));
    end

    function val = getGradientLinear(xj, t_now, gForm, ...
                p_max, p_min, gStart, gEnd, v, r_off, sWidth, b2g)
        if t_now < b2g, t_eff = 0; else, t_eff = t_now - b2g; end
        switch gForm
            case 1
                if xj <= gStart,       val = p_min;
                elseif xj >= gEnd,     val = p_max;
                else,  val = p_min + (p_max-p_min)*(xj-gStart)/(gEnd-gStart);
                end
            case {2, 3}
                step_center = gStart + v * t_eff;
                left  = step_center - sWidth/2;
                right = step_center + sWidth/2;
                if xj < left,          val = p_min;
                elseif xj > right,     val = p_max;
                else,  val = p_min + (p_max-p_min)*(xj-left)/sWidth;
                end
            otherwise,                 val = p_max;
        end
    end

    function val = getGradientConcave(xj, t_now, gForm, ...
                p_max, p_min, gStart, gEnd, v, r_off, p, sWidth, t_frac, y_thr, b2g)
        if t_now < b2g, t_eff = 0; else, t_eff = t_now - b2g; end
        switch gForm
            case 1
                if xj <= gStart,       val = p_min;
                elseif xj >= gEnd,     val = p_max;
                else
                    x_norm = (xj-gStart)/(gEnd-gStart);
                    val    = p_min + (p_max-p_min)*x_norm^p;
                end
            case 2
                step_center = r_off + v * t_eff;
                left  = step_center - sWidth/2;
                right = step_center + sWidth/2;
                if xj < left,          val = p_min;
                elseif xj > right,     val = p_max;
                else
                    x_norm = (xj-left)/sWidth;
                    val    = p_min + (p_max-p_min)*x_norm^p;
                end
            case 3
                step_center = gStart + v * t_eff;
                left  = step_center - sWidth/2;
                x_cut = left + y_thr * sWidth;
                if xj < x_cut,         val = p_min;
                else,                  val = p_max;
                end
            otherwise,                 val = p_max;
        end
    end

end