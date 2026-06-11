
% main.m
% Find ground state and run a simple time evolution (split-step FFT).
clearvars; close all; clc;

% ---------------- User Parameters ----------------
rspan = [1e-6, 10];                  % radial domain [r0, R]
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@asymptotic_event);

nth_state = 1;                       % 1 = ground state
Q_left_init = 0;
Q_right_init = 15;
tol_target = 1e-8;

% Time-evolution / spectral grid parameters
Mx = 256; My = 256;                  % grid resolution
a  = 20;                             % domain [-a,a] x [-a,a]
dt = 0.01;                           % time step
T  = 10;                            % final time
epsilon = -1;                        % NLSE nonlinearity sign

doTimeEvolution = true;              % set false to skip evolution

% Output folder
results_dir = 'results';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

% Optional: adjust event tolerances (useful for very tight runs)
% NLSE_EVENT_CFG.epsQ = 1e-8;
% NLSE_EVENT_CFG.epsP = 1e-8;
% assignin('base','NLSE_EVENT_CFG', NLSE_EVENT_CFG);

% ---------------- Find Ground State ----------------
fprintf('Finding ground state (nth_state = %d)...\n', nth_state);
[r_sol, Q_sol, Q_star, diagnostics] = find_nth_state(nth_state, opts, rspan, Q_left_init, Q_right_init, tol_target);
fprintf('Found Q(0) ≈ %.12g\n', Q_star);

% Save radial solution
save(fullfile(results_dir, sprintf('radial_n%d_%s.mat', nth_state, timestamp)), 'r_sol', 'Q_sol', 'Q_star', 'diagnostics', '-v7.3');

% Plot radial solution
figure('Name','Radial Solution','NumberTitle','off');
plot(r_sol, Q_sol, 'LineWidth', 2);
xlabel('r'); ylabel('Q(r)'); grid on;
title(sprintf('Ground State Q(0)=%.6g', Q_star));
saveas(gcf, fullfile(results_dir, sprintf('radial_plot_n%d_%s.png', nth_state, timestamp)));

% ---------------- Build 2D Initial Condition ----------------
dx = 2*a / Mx;
dy = 2*a / My;
x = linspace(-a, a-dx, Mx);
y = linspace(-a, a-dy, My);
[X, Y] = meshgrid(x, y);
R_grid = sqrt(X.^2 + Y.^2);

% Interpolate radial solution onto 2D grid (spline; zero outside R)
Q_on_R = interp1(r_sol, Q_sol, R_grid, 'spline', 0);
psi0 = Q_on_R;   % real-valued initial field

% Visualize initial |psi|^2
figure('Name','Initial |psi|^2','NumberTitle','off');
h1 = surf(X, Y, abs(psi0).^2, 'EdgeColor', 'none');
shading interp; colorbar;
xlabel('x'); ylabel('y'); zlabel('|psi|^2');
title('Initial Condition |psi|^2');
xlim([-5,5]); ylim([-5,5]);
saveas(gcf, fullfile(results_dir, sprintf('psi0_surface_n%d_%s.png', nth_state, timestamp)));

% Save psi0
save(fullfile(results_dir, sprintf('psi0_n%d_%s.mat', nth_state, timestamp)), 'psi0', 'X', 'Y', '-v7.3');

% ---------------- Time Evolution (with saving for replay) ----------------
if doTimeEvolution
    cfg = struct();
    cfg.Mx = Mx; cfg.My = My; cfg.a = a; cfg.dt = dt; cfg.T = T; cfg.epsilon = epsilon;
    fprintf('Running time evolution to T = %.4g ...\n', T);

    % Video and frame save settings
    video_fname = fullfile(results_dir, sprintf('timeevo_movie_n%d_%s.mp4', nth_state, timestamp));
    mat_fname   = fullfile(results_dir, sprintf('timeevo_frames_n%d_%s.mat', nth_state, timestamp));
    fps = 10;                % frames per second in output video
    Nskip = max(1, round(1/(dt*fps))); % how many time steps per saved frame
    maxSavedFrames = 500;    % protect against huge .mat files

    % Initialize video writer
    v = VideoWriter(video_fname, 'MPEG-4'); % or 'Motion JPEG AVI' if MP4 not supported
    v.FrameRate = fps;
    open(v);

    % Preallocate small frame cache (will expand if needed). Store grayscale frames.
    frames = zeros(min(maxSavedFrames, ceil(T/dt/Nskip)), Mx, My, 'uint8');
    saved_idx = 0;

    % Prepare spectral operator
    mx = (-Mx/2 : Mx/2-1);
    my = (-My/2 : My/2-1);
    kx = (pi/a) * mx;
    ky = (pi/a) * my;
    [KX, KY] = meshgrid(kx, ky);
    lambda = KX.^2 + KY.^2;
    K = exp(-1i * lambda * dt/2);

    % Initialize psi from radial Q_sol (already computed above)
    psi = Q_on_R;

    nt = ceil(T / dt);
    t = 0;
    frameCount = 0;

    % For progress display
    tic;
    for it = 1:nt
        % Strang split-step
        psi_hat = fftshift(fft2(psi));
        psi_hat = K .* psi_hat;
        psi = ifft2(ifftshift(psi_hat));

        psi = exp(-1i * epsilon * dt * abs(psi).^2) .* psi;

        psi_hat = fftshift(fft2(psi));
        psi_hat = K .* psi_hat;
        psi = ifft2(ifftshift(psi_hat));

        t = t + dt;

        % Save a frame every Nskip steps
        if mod(it, Nskip) == 0
            frameCount = frameCount + 1;

            % Create an RGB image from |psi|^2 using colormap
            I = abs(psi).^2;
            I = I - min(I(:)); I = I / max(I(:) + eps); % normalize to [0,1]
            cmap = hot(256);
            RGB = ind2rgb( round(I*(size(cmap,1)-1)), cmap );

            % Write video frame
            writeVideo(v, im2frame(RGB));

            % Save to frames cache (grayscale uint8) if under limit
            if saved_idx < size(frames,1)
                saved_idx = saved_idx + 1;
                % store resized/grayscale to reduce size if needed
                G = im2uint8(mat2gray(I));
                frames(saved_idx,:,:) = G;
            end
        end

        % Optional: quick console progress
        if mod(it, max(1,round(nt/10))) == 0
            fprintf('Progress: %.0f%% (t = %.3g / %.3g)\n', 100*it/nt, t, T);
        end
    end
    elapsed = toc;
    fprintf('Time evolution finished in %.2f s, %d frames written.\n', elapsed, frameCount);

    close(v); % finalize video file

    % Trim frames cache to actual saved size
    frames = frames(1:saved_idx,:,:);

    % Save final field and diagnostics + frames
    phase_info = struct(); % compute lightweight phase_info if desired
    % Quick recomputation of phase and mass over movie frames (optional)
    % Here we save simple placeholders; you can compute full phase_info inside loop if needed.
    save(mat_fname, 'psi', 'frames', 'phase_info', '-v7.3');

    % Save final field and diagnostics (as before)
    save(fullfile(results_dir, sprintf('timeevo_n%d_%s.mat', nth_state, timestamp)), 'psi', 'phase_info', '-v7.3');

    % Visualize final |psi|^2
    figure('Name','Final |psi|^2','NumberTitle','off');
    surf(X, Y, abs(psi).^2, 'EdgeColor', 'none');
    shading interp; colorbar;
    xlabel('x'); ylabel('y'); zlabel('|psi|^2');
    title(sprintf('Final |psi|^2 at T=%.3g', T));
    xlim([-5,5]); ylim([-5,5]);
    saveas(gcf, fullfile(results_dir, sprintf('psi_final_surface_n%d_%s.png', nth_state, timestamp)));

    % Optionally also save a small preview image
    imwrite(RGB, fullfile(results_dir, sprintf('preview_final_n%d_%s.png', nth_state, timestamp)));

    fprintf('Saved video: %s\nSaved frames MAT: %s\n', video_fname, mat_fname);
end
fprintf('Main run complete. Results saved to: %s\n', results_dir);