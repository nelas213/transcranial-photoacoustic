
% =========================================================================
% SIMULATION
% =========================================================================

load('C:\Users\tpas\Desktop\Dateidatei\Simple Segmentation\outputs simpa\final_times_series_data_0.25_pitch.mat')

% change scale to 2 to reproduce the higher resolution figures used in the
% help file
scale = 1;

input_file=input_file
output_file=output_file

% create the computational grid                  % size of the PML in grid points
Nx = 300;
Ny = 300;
Nz = 410;
dx = 0.125e-3;            % grid point spacing in the x direction [m]
dy = 0.125e-3;            % grid point spacing in the y direction [m]
dz = 0.125e-3;            % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
PML_size = getOptimalPMLSize([Nx,Ny,Nz]); 

% define the properties of the propagation medium
load('input_file/sos_perfect_depth_150.mat')
medium.sound_speed = sos;      % [m/s]
%load('input_file/sos_perfect_depth_400.mat')
%medium.sound_speed = sos_tief;
medium.density = 1000*ones(Nx, Ny, Nz);
medium.alpha_coeff = 0.01;

medium.alpha_power = 1.1;
medium.alpha_mode = 'no_dispersion';

%%% 
% define a binary planar sensor
sensor.mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
sensor.mask(1:2:300, 1:2:300, 1) = 1;

% create the time array
%kgrid.setTime(1612,1/4e7);
% load sampling rate from settings
dt = 1.0 / double(40 * 1000000);

% Simulate as many time steps as a wave takes to traverse diagonally through the entire tissue
Nt = round((sqrt(Ny*Ny+Nx*Nx+Nz*Nz)*dx / mean(medium.sound_speed, 'all')) / dt);

estimated_cfl_number = dt / dx * mean(medium.sound_speed, 'all');

% smaller time steps are better for numerical stability in time progressing simulations
% A minimum CFL of 0.3 is advised in the kwave handbook.
% In case we specify something larger, we use a higher sampling rate than anticipated.
% Otherwise we simulate with the target sampling rate
if estimated_cfl_number < 0.3
    kgrid.setTime(Nt, dt);
else
    kgrid.t_array = makeTime(kgrid, medium.sound_speed, 0.3);
end

% % load sampling rate from settings
% dt = 1.0 / double(40 * 1000000);
% % 
% % % Simulate as many time steps as a wave takes to traverse diagonally through the entire tissue
% Nt = 3622;
% % 
% kgrid.setTime(Nt, 13.5365/1000000000);




% set the input arguements
input_args = {'PMLSize', PML_size, 'PMLInside', false, ...
    'PlotPML', false, 'Smooth', false, 'DataCast', 'single'};

% %load input data as sensor data p0 only vessel
load('input_file/depth_150_larger_p0_just_blood.mat',"sensor_data")
%load('input_file/depth_400_larger_p0_just_blood.mat',"sensor_data")

% load data for the whole p0
% %load('input_file/depth_150_larger_p0.mat',"sensor_data")
% load('input_file/depth_400_larger_p0.mat',"sensor_data")


% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

source.p0=0;

% run the time-reversal reconstruction
p0_recon = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});
p0_recon = gather(p0_recon);

% apply a positivity condition
% p0_recon(p0_recon < 0) = 0;

% =========================================================================
% VISUALISATION
% =========================================================================

%% 
% plot the reconstructed initial pressure
figure;
subplot(2, 2, 1);
imagesc(max(p0_recon(:, :, :),  [], 3));
title('x-y plane');
axis image;

subplot(2, 2, 2);
imagesc(squeeze(p0_recon(:, 92, :)));
title('x-z plane');
axis image;
xlabel('(All axes in mm)');

subplot(2, 2, 3);
imagesc(squeeze(p0_recon(92, :, :)));
title('y-z plane');
axis image;
colormap(getColorMap);

%% 

% ====================================================
%   SAVE
% =====================================================
% save for the p0 just blood
save('output_file/perfect_sos_depth_150_larger_p0_just_blood.mat',"p0_recon")
%save('output_file/perfect_sos_depth_400_larger_p0_just_blood.mat',"p0_recon")

% save for the whole p0
% %save('output_file/perfect_sos_depth_150_larger_p0.mat',"p0_recon")
% save('output_file/perfect_sos_depth_400_larger_p0.mat',"p0_recon")