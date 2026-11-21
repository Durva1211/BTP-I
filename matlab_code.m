clear; clc;

%% 1. DESIGN (from the Python output)
name      = "singlet_doublet N-LASF9/N-FK51A/N-SF6";
EFL       = 40.0;
fnum      = 3.6;
psf_c     = 1.02;
psf_e     = 0.8;
psf_worst = 0.99;
L_optics  = 43.20;
tau       = 0.932;

%% 2. SENSOR (HAS3)
N_pix = 1280; p_mm = 0.011;
QE = 0.60; read_noise = 5.0; dark_rate = 5.0; sky_rate = 2.0;   % e-, e-/s, e-/s

%% 3. REQUIREMENTS
m_lim = 6.0; t_int = 0.050;
SNR_det = 5; SNR_use = 10; min_stars = 6;
max_cross = 5.0;                    % [arcsec]
sun_excl = 45;                      % baffle exclusion angle [deg]
sig_thermal = 1.2; sig_calib = 1.2; % [arcsec]
sig_centroid_min = 0.08;            % best centroid accuracy achievable [px]

%% 4. PHOTOMETRY AND 1U BUDGET
F0 = 3e10;                          % photons/s/m^2 from an Mv 0 star (sensor band)
density_derate = 0.5;               % worst-case star density
U = 100; pcb = 30; fasteners = 6; vane_margin = 3; min_spare = 3;

%% ---------------- calculations ----------------
D_ap    = EFL / fnum;                              % aperture [mm]
pix_as  = p_mm / EFL * 206264.8;                   % arcsec per pixel
FOV     = 2 * atand(N_pix * p_mm / 2 / EFL);       % square FOV [deg]
FOV_d   = 2 * atand(sqrt(2) * N_pix * p_mm / 2 / EFL);   % diagonal FOV [deg]

% light from the faintest star [e-]
A = pi * (D_ap * 1e-3 / 2)^2;
N_star_e = F0 * 10^(-0.4 * m_lim) * A * tau * QE * t_int;
noise_bg = read_noise^2 + (dark_rate + sky_rate) * t_int;

% peak pixel: average of star centred on a pixel and on a pixel corner
peak_frac = @(s) 0.5 * (erf(0.5 ./ (sqrt(2)*s)).^2 + (0.5*erf(1 ./ (sqrt(2)*s))).^2);
snr = @(s) N_star_e * peak_frac(s) ./ sqrt(N_star_e * peak_frac(s) + noise_bg);

% faintest magnitude reaching a given SNR
peak_needed = @(k) (k^2 + sqrt(k^4 + 4*k^2*noise_bg)) / 2;
mag_for = @(s, k) m_lim - 2.5 * log10(peak_needed(k) / (N_star_e * peak_frac(s)));

% stars in view (all-sky count brighter than m: 6.57 e^(1.08 m))
n_stars = @(m) density_derate * 6.57 * exp(1.08 * m) / 41253 * FOV^2;

% baffle and 1U fit
baffle = D_ap / (tand(sun_excl) - tand(FOV_d / 2)) + vane_margin;
spare  = U - pcb - fasteners - L_optics - baffle;

%% ---------------- report ----------------
fprintf('\n=== %s | EFL %.1f mm f/%.1f ===\n', name, EFL, fnum);
fprintf('Aperture %.2f mm | %.1f arcsec/px | FOV %.1f deg (diag %.1f)\n', D_ap, pix_as, FOV, FOV_d);
fprintf('Signal from Mv %.0f star: %.0f e-\n', m_lim, N_star_e);
fprintf('\nLength: optics %.1f + baffle %.1f + PCB %d + fasteners %d = %.1f mm -> spare %.1f mm\n', ...
    L_optics, baffle, pcb, fasteners, U - spare, spare);

cases = ["Centre", "Edge", "Worst"];
sigs  = [psf_c, psf_e, psf_worst];
fprintf('\n%-7s %6s %8s %8s %7s %9s\n', 'PSF', 'sigma', 'SNR@Mv6', 'Mag lim', 'Stars', 'Cross["]');
pass_all = true;
for i = 1:3
    s = sigs(i);
    m_det = mag_for(s, SNR_det);
    m_use = mag_for(s, SNR_use);
    n_det = n_stars(m_det);
    n_use = max(n_stars(m_use), 1);
    sig_c = max(sig_centroid_min, s / SNR_use);          % centroid error [px]
    cross = sqrt((pix_as * sig_c / sqrt(n_use))^2 + sig_thermal^2 + sig_calib^2);
    ok = snr(s) >= SNR_det && n_det >= min_stars && cross <= max_cross;
    pass_all = pass_all && ok;
    fprintf('%-7s %6.2f %8.2f %8.2f %7.1f %9.2f  %s\n', cases(i), s, snr(s), m_det, n_det, cross, ...
        pick(ok, 'ok', 'FAIL'));
end

fits = spare >= min_spare;
fprintf('\nOptics: %s | 1U fit: %s\n', pick(pass_all, 'PASS', 'FAIL'), pick(fits, 'PASS', 'FAIL'));
if pass_all && fits
    fprintf('DESIGN PASSES\n');
else
    fprintf('DESIGN FAILS\n');
end

function out = pick(b, a, c)
    if b, out = a; else, out = c; end
end