function sweeping_pulse_spacetime()
%% sweeping_pulse_spacetime.m
%
% Simulates a Gaussian transcriptional pulse (intronic signal) travelling
% at constant speed from posterior to anterior along the AP axis.
% Records intronic and exonic signals at a fixed posterior observation
% point to show how a spatial wave produces a temporal pulse.
%
% Outputs:
%   - Animated movie:  sweeping_pulse_spatial_profiles.mp4
%   - PDF Figure 1:    sweeping_pulse_spatial_profiles.pdf
%   - PDF Figure 2:    sweeping_pulse_posterior_timecourse.pdf
%
% Color scheme:
%   Intronic  -->  gray  [0.5 0.5 0.5]
%   Exonic    -->  black [0   0   0  ]

%% ── Parameters ───────────────────────────────────────────────────────────────
L      = 1;
Nx     = 300;
z      = linspace(0, L, Nx)';
T      = 5;
dt     = 0.01;
tvec   = 0:dt:T;
Nt     = numel(tvec);
lambda = 1.0;
A      = 1.0;
sigma  = 0.08;
v      = 0.2;          % constant sweep speed

%% ── Color scheme ─────────────────────────────────────────────────────────────
col_intronic = [0.5 0.5 0.5];   % gray
col_exonic   = [0   0   0  ];   % black

%% ── Pulse starts outside the domain ─────────────────────────────────────────
z_start = L + 3*sigma;

%% ── Observation point (posterior end) ───────────────────────────────────────
z_obs = 1;
[~, iz] = min(abs(z - z_obs));

%% ── Initialise signals ───────────────────────────────────────────────────────
x    = zeros(Nx, 1);   % intronic
y    = zeros(Nx, 1);   % exonic
x_tc = zeros(Nt, 1);   % intronic time course at z_obs
y_tc = zeros(Nt, 1);   % exonic  time course at z_obs

%% ── Figure 1: Animated spatial profiles ─────────────────────────────────────
fig1 = figure('Color', 'w', 'Name', 'Sweeping Pulse — Spatial Profiles');

hInt = plot(z, zeros(Nx,1), '-',  'Color', col_intronic, 'LineWidth', 2);
hold on;
hExo = plot(z, zeros(Nx,1), '-',  'Color', col_exonic,   'LineWidth', 2);
hold off;

axis([0 L 0 1.2*A]);
xlabel('AP position (posterior \rightarrow anterior)');
ylabel('Expression');
legend({'intronic', 'exonic (norm)'}, 'Location', 'northwest');
legend('boxoff');
set(gca, 'FontSize', 12);
grid on; box off;

%% ── Movie writer setup ───────────────────────────────────────────────────────
movie_file   = 'sweeping_pulse_spatial_profiles.mp4';
frame_interval = 5;
vw = VideoWriter(movie_file, 'MPEG-4');
vw.FrameRate = max(1, round(1 / (dt * frame_interval)));
open(vw);

%% ── Time integration (Euler) + movie capture ─────────────────────────────────
for k = 1:Nt
    t = tvec(k);

    % Pulse centre moves from posterior to anterior
    z0 = z_start - v*t;

    % Intronic signal (prescribed Gaussian pulse)
    x = A * exp(-(z - z0).^2 / (2*sigma^2));

    % Exonic dynamics (first-order ODE, Euler step)
    y = y + dt*(x - lambda*y);

    % Record time courses at observation point
    x_tc(k) = x(iz);
    y_tc(k) = y(iz);

    % Visual normalisation of exonic for spatial plot
    max_y = max(y);
    if max_y > 0
        y_plot = y * (max(x) / max_y);
    else
        y_plot = y;
    end

    % Update animated plot
    set(hInt, 'YData', x);
    set(hExo, 'YData', y_plot);
    title(sprintf('Sweeping pulse — spatial profiles,  t = %.2f', t));
    drawnow;

    % Capture frame
    if mod(k, frame_interval) == 0
        writeVideo(vw, getframe(fig1));
    end
end

close(vw);
fprintf('Movie saved: %s\n', movie_file);

% Save Figure 1 as PDF (final frame)
exportgraphics(fig1, 'sweeping_pulse_spatial_profiles.pdf', 'ContentType', 'vector');
fprintf('Figure saved: sweeping_pulse_spatial_profiles.pdf\n');

%% ── Global normalisation of exonic time course ───────────────────────────────
max_x_tc = max(x_tc);
max_y_tc = max(y_tc);
if max_y_tc > 0
    scale_tc = max_x_tc / max_y_tc;
else
    scale_tc = 1;
end
y_tc_norm = y_tc * scale_tc;

%% ── Figure 2: Posterior time courses ────────────────────────────────────────
fig2 = figure('Color', 'w', 'Name', 'Sweeping Pulse — Posterior Time Course');

plot(tvec, x_tc,      '-', 'Color', col_intronic, 'LineWidth', 2); hold on;
plot(tvec, y_tc_norm, '-', 'Color', col_exonic,   'LineWidth', 2);
hold off;

xlabel('Time');
ylabel('Normalised expression');
title(sprintf('Posterior time course at z = %.2f', z_obs));
legend({'intronic', 'exonic (norm)'}, 'Location', 'northwest');
legend('boxoff');
grid on; box off;
set(gca, 'FontSize', 12);

% Save Figure 2 as PDF
exportgraphics(fig2, 'sweeping_pulse_posterior_timecourse.pdf', 'ContentType', 'vector');
fprintf('Figure saved: sweeping_pulse_posterior_timecourse.pdf\n');

end