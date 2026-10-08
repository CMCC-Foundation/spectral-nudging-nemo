% =========================================================================
% generate_gamma_map.m
%
% Generates the 3D spatially-varying nudging coefficient (Gamma map) for
% a regional model for the Adriatic Sea (1/16 degree, NEMO v3.6).
%
% This script is provided as an example of how the Gamma map was generated
% for the Adriatic configuration. It is domain-specific and will need to
% be adapted for other configurations.
%
% The Gamma map defines where and how strongly the spectral nudging is
% applied. It ranges from 0 (no nudging) to 1 (maximum nudging) and is
% defined as the product of:
%   - a horizontal 2D field based on bathymetry:
%       deeper ocean = stronger nudging (via hyperbolic tangent function)
%   - a vertical 1D profile:
%       nudging is zero near the surface, maximum at a reference depth,
%       and zero at the bottom
%
% The resulting 3D field is saved as Gamma_map_3d.nc, which is specified
% in the namelist as cn_gamma_file.
%
% NOTE: this script requires the x_y_remove.mat file (provided in the
% same directory) which contains manually selected grid points used to
% zero out the Gamma map along the eastern Adriatic coastline. This was
% done to avoid nudging at specific coastal areas.
%
% =========================================================================

clear all; close all; clc;

% =========================================================================
%% USER SETTINGS - adapt paths to your local directory
% =========================================================================

% Input files
bath_file = './depth_ADRb.060_04.nc';   % NEMO bathymetry file
mask_file = './mask.nc';             % NEMO mask file

% Variable names
bath_var  = 'Bathymetry';
tmask_var = 'tmask';
depth_var = 'gdept_0';
lon_var   = 'nav_lon';
lat_var   = 'nav_lat';

% Output file
output_file = 'Gamma_map_3d.nc';

% =========================================================================
%% HORIZONTAL GAMMA PARAMETERS
% These were tuned for the Adriatic configuration.
% Adapt for your domain based on your bathymetry and scientific objectives.
% =========================================================================

ho       = 400;  % (m) reference depth for hyperbolic tangent transition
                 % areas deeper than ho will have Gamma closer to 1
                 % areas shallower than ho will have Gamma closer to 0
dh       = 300;  % (m) width of the transition zone around ho
                 % larger dh = smoother transition between nudged/non-nudged areas
min_bath = 40;   % (m) minimum bathymetry for nudging
                 % grid points shallower than this will have Gamma = 0

% =========================================================================
%% VERTICAL GAMMA PARAMETERS
% These were tuned for the AdriaCLIM configuration.
% =========================================================================

min_surf  = 5;    % (m) depth above which Gamma is zero
                  % no nudging applied near the surface
max_depth = 150;  % (m) depth of maximum Gamma (= 1)
                  % Gamma increases from 0 at min_surf to 1 at max_depth
                  % then decreases back to 0 at the bottom

% =========================================================================
%% READ INPUT FILES
% =========================================================================

bath  = double(ncread(bath_file, bath_var));
lon   = ncread(bath_file, lon_var);
lat   = ncread(bath_file, lat_var);

tmask = double(ncread(mask_file, tmask_var));
dz    = ncread(mask_file, depth_var);
dz1d  = squeeze(dz(1,1,:));   % 1D depth levels

nx = size(bath, 1);
ny = size(bath, 2);
nz = size(tmask, 3);

fprintf('Grid dimensions: %d x %d x %d\n', nx, ny, nz);

% =========================================================================
%% HORIZONTAL GAMMA (2D)
% =========================================================================

%--- Plot raw bathymetry for reference
bath_plot = bath; bath_plot(bath == 0) = NaN;
figure(1);
pcolor(bath_plot'); shading flat; colorbar; colormap('jet');
title('Bathymetry (m)'); xlabel('i'); ylabel('j');

%--- Smooth bathymetry with a 5x5 box average
%    to avoid sharp transitions at steep topography
bath_smooth = conv2(bath, ones(5)/5^2, 'same');

figure(2);
bath_smooth_plot = bath_smooth; bath_smooth_plot(bath_smooth == 0) = NaN;
pcolor(bath_smooth_plot'); shading flat; colorbar; colormap('jet');
title('Smoothed bathymetry (m)'); xlabel('i'); ylabel('j');

%--- Apply hyperbolic tangent to define horizontal Gamma
%    Go2d ranges from ~0 (shallow) to ~1 (deep)
%    transition centered at ho with width controlled by dh
Go2d = (1 + tanh((bath_smooth - ho) ./ dh)) ./ 2;

%--- Remove very shallow areas
Go2d(bath < min_bath) = 0;

%--- Apply 2D ocean mask (land = 0)
mask2d = squeeze(tmask(:,:,1));
Go2d   = Go2d .* mask2d;

%--- Fix boundary columns to avoid edge effects
%    replicate the third column into the first two
for i = 1:2
    Go2d(:,i) = Go2d(:,3);
end

figure(3);
pcolor(Go2d'); shading flat; colorbar; colormap('jet'); caxis([0 1]);
title('Horizontal Gamma (2D) before boundary adjustments');
xlabel('i'); ylabel('j');

%--- Adriatic-specific: zero out eastern coastline features
%    x_y_remove.mat contains manually selected grid points (using ginput)
%    that define the boundary above/east of which Gamma is set to zero.
%    This was done to avoid nudging in shallow coastal areas.
load('x_y_remove.mat')
for ii = 1:length(x)
    Go2d(round(x(ii):end), round(y(ii):end)) = 0;
end

figure(4);
pcolor(Go2d'); shading flat; colorbar; colormap('jet'); caxis([0 1]);
title('Horizontal Gamma (2D) after boundary adjustments');
xlabel('i'); ylabel('j');

% =========================================================================
%% VERTICAL GAMMA PROFILE (1D)
% =========================================================================

% Find depth level indices closest to user-defined depths
[~, i_max]       = min(abs(dz1d - max_depth));
[~, i_min_surf]  = min(abs(dz1d - min_surf));
[~, i_min_depth] = min(abs(dz1d - dz1d(end)));

fprintf('\nVertical Gamma profile:\n');
fprintf('   Zero above %.1f m (level %d)\n', dz1d(i_min_surf), i_min_surf);
fprintf('   Maximum (=1) at %.1f m (level %d)\n', dz1d(i_max), i_max);
fprintf('   Zero at bottom %.1f m (level %d)\n', dz1d(i_min_depth), i_min_depth);

% Define anchor points for vertical Gamma profile
gamma_1d               = NaN(nz, 1);
gamma_1d(i_max)        = 1;   % maximum nudging at max_depth
gamma_1d(1:i_min_surf) = 0;   % no nudging near surface
gamma_1d(i_min_depth)  = 0;   % no nudging at bottom

% Interpolate between anchor points using pchip
gamma_1d = fillmissing(gamma_1d, 'pchip');

figure(5);
plot(gamma_1d, dz1d, 'b-o', 'LineWidth', 2, 'MarkerSize', 4);
set(gca, 'YDir', 'reverse');
xlabel('Gamma'); ylabel('Depth (m)');
title('Vertical Gamma profile');
grid on;
xlim([0 1.1]);

% =========================================================================
%% COMBINE HORIZONTAL AND VERTICAL GAMMA (3D)
% =========================================================================

Go3d = zeros(nx, ny, nz);
for k = 1:nz
    Go3d(:,:,k) = Go2d .* gamma_1d(k);
end

% Apply 3D ocean mask
Go3d = Go3d .* tmask;

%--- Plot 3D Gamma at selected depth levels
depth_plot = [1, 20, 40, 80];   % level indices to plot
figure(6);
for p = 1:length(depth_plot)
    subplot(2, 2, p);
    pcolor(Go3d(:,:,depth_plot(p))'); shading flat;
    colorbar; colormap('jet'); caxis([0 1]);
    title(sprintf('Gamma at %.0f m', dz1d(depth_plot(p))));
    xlabel('i'); ylabel('j');
end
sgtitle('3D Gamma map at selected depth levels');

% =========================================================================
%% SAVE OUTPUT
% =========================================================================

if exist(output_file, 'file')
    delete(output_file);
end

nccreate(output_file, 'Go', ...
    'Dimensions', {'x', nx, 'y', ny, 'z', nz}, ...
    'Datatype', 'double', ...
    'Format', 'classic');
ncwrite(output_file, 'Go', Go3d);

% Add variable attributes
ncwriteatt(output_file, 'Go', 'long_name', 'Spectral nudging coefficient (Gamma map)');
ncwriteatt(output_file, 'Go', 'units', '1');
ncwriteatt(output_file, 'Go', 'valid_min', 0.0);
ncwriteatt(output_file, 'Go', 'valid_max', 1.0);

% Add global attributes documenting the parameters used
ncwriteatt(output_file, '/', 'history',      ['Created by generate_gamma_map.m on ' datestr(now)]);
ncwriteatt(output_file, '/', 'configuration','Adriatic Sea 1/16 degree NEMO v3.6');
ncwriteatt(output_file, '/', 'ho_m',         ho);
ncwriteatt(output_file, '/', 'dh_m',         dh);
ncwriteatt(output_file, '/', 'min_bath_m',   min_bath);
ncwriteatt(output_file, '/', 'min_surf_m',   min_surf);
ncwriteatt(output_file, '/', 'max_depth_m',  max_depth);

fprintf('\nGamma map saved to %s\n', output_file);
fprintf('Variable: Go, dimensions: %d x %d x %d\n', nx, ny, nz);
