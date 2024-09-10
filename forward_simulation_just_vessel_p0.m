% =========================================================================
% SIMULATION
% =========================================================================

% change scale to 2 to reproduce the higher resolution figures used in the
% help file
scale = 1;

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
%sos
load('C:/Users/tpas/Desktop/Dateidatei/final work nelas/medium/sos_perfect_depth_150.mat')
%load('C:/Users/tpas/Desktop/Dateidatei/final work nelas/medium/sos_perfect_depth_400.mat')
medium.sound_speed = sos;
% density
load('C:/Users/tpas/Desktop/Dateidatei/final work nelas/medium/density_perfect_depth_150.mat')
%load('C:/Users/tpas/Desktop/Dateidatei/final work nelas/medium/density_perfect_depth_400.mat')
medium.density = density;
%attenuation
load('C:/Users/tpas/Desktop/Dateidatei/final work nelas/medium/alpha_0_perfect_depth_150.mat')
%load('C:/Users/tpas/Desktop/Dateidatei/final work nelas/medium/alpha_0_perfect_depth_400.mat')
medium.alpha_coeff = alpha_0;

medium.alpha_power = 1.1;
medium.alpha_mode = 'no_dispersion';


% assign to the source structure geringe tiefe
load('C:/Users/tpas/Desktop/Dateidatei/final work nelas/mat file p0/p0_150_just_blood.mat')
source.p0 = p0_150_inPa_blood;

% assign to the source structure hohe tiefe
%load('C:/Users/tpas/Desktop/Dateidatei/mat file p0/p0_400.mat')
%source.p0 = p0_400_inPa;


% define a binary planar sensor
sensor.mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz); 
sensor.mask(1:2:300, 1:2:300, 1) = 1;

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

% set the input arguements
input_args = {'PMLSize', PML_size, 'PMLInside', false, ...
    'PlotPML', false, 'Smooth', false, 'DataCast', 'single'};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
sensor_data = gather(sensor_data);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          save 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('C:/Users/tpas/Desktop/Dateidatei/final work nelas/time series data/larger p0/depth_150_larger_p0_just_blood.mat',"sensor_data")


clear all

% % =========================================================================
% % SIMULATION FOR DEPTH 400
% % =========================================================================
% 
% =========================================================================
% SIMULATION
% =========================================================================

% change scale to 2 to reproduce the higher resolution figures used in the
% help file
scale = 1;

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
%sos
%load('C:/Users/tpas/Desktop/Dateidatei/final work nelas/medium/sos_perfect_depth_150.mat')
load('C:/Users/tpas/Desktop/Dateidatei/final work nelas/medium/sos_perfect_depth_400.mat')
medium.sound_speed = sos_tief;
% density
%load('C:/Users/tpas/Desktop/Dateidatei/final work nelas/medium/density_perfect_depth_150.mat')
load('C:/Users/tpas/Desktop/Dateidatei/final work nelas/medium/density_perfect_depth_400.mat')
medium.density = density_tief;
%attenuation
%load('C:/Users/tpas/Desktop/Dateidatei/final work nelas/medium/alpha_0_perfect_depth_150.mat')
load('C:/Users/tpas/Desktop/Dateidatei/final work nelas/medium/alpha_0_perfect_depth_400.mat')
medium.alpha_coeff = alpha_0_tief;

medium.alpha_power = 1.1;
medium.alpha_mode = 'no_dispersion';


% assign to the source structure geringe tiefe
load('C:/Users/tpas/Desktop/Dateidatei/final work nelas/mat file p0/p0_400_just_blood.mat')
source.p0 = p0_400_inPa_blood;

% assign to the source structure hohe tiefe
%load('C:/Users/tpas/Desktop/Dateidatei/mat file p0/p0_400.mat')
%source.p0 = p0_400_inPa;

%================================================
%Hier habe ich zu wenig Detektoren kann da eigentlich bis 300 gehen
%=================================================================



% define a binary planar sensor
sensor.mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
sensor.mask(1:2:300, 1:2:300, 1) = 1;

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

% set the input arguements
input_args = {'PMLSize', PML_size, 'PMLInside', false, ...
    'PlotPML', false, 'Smooth', false, 'DataCast', 'single'};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
sensor_data = gather(sensor_data);
% 
save('C:/Users/tpas/Desktop/Dateidatei/final work nelas/time series data/larger p0/depth_400_larger_p0_just_blood.mat',"sensor_data")