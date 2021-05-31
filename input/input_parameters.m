%==========================================================================
% set basic simulation parameters
%==========================================================================

Lx=25e6;     % model extension in x-direction [m]
Lz=15e6;     % model extension in y-direction [m]

nx=360;     % grid points in x-direction
nz=260;     % grid points in z-direction

dt=5.;     % time step [s]
nt=1000;     % number of iterations

order=4;    % finite-difference order (2 or 4)

%==========================================================================
% model type
%==========================================================================

model_type=1;

% 1=homogeneous 
% 2=homogeneous with localised density perturbation
% 3=layered medium
% 4=layered with localised density perturbation
% 5=vertical gradient medium
% 6=vertical gradient medium with localised density perturbation

% 'initial'= read initial model for waveform inversion (mu_initial, rho_initial)

%==========================================================================
% source-time function
%==========================================================================

f_min=0.005;     % minimum frequency [Hz]
f_max=0.01;     % maximum frequency [Hz]

%==========================================================================
% simulation mode
%==========================================================================

simulation_mode='noise_source_kernels';

% 'forward'                 regular forward simulation
% 'forward_green'           forward simulation where Fourier transform of the Greens function is computed on-the-fly (preparation to compute correlation function)
% 'correlation'             compute correlation functions as output, requires Fourier-transformed Green function to be present
% 'noise_source_kernel'     compute sensitivity kernel for the noise source power-spectral density distribution

%==========================================================================
% source positions
%==========================================================================

%src_x=[15.83341e6];
src_x=[13.5e6];
src_z=[7.5e6];

%==========================================================================
% receiver positions
%==========================================================================

%- just one receiver
%rec_x=[9.16659e6];
rec_x=[11.5e6];
rec_z=[7.5e6];

%==========================================================================
% absorbing boundaries
%==========================================================================

width=2e6;     % width of the boundary layer in m

absorb_left=1;  % absorb waves on the left boundary
absorb_right=1; % absorb waves on the right boundary
absorb_top=1;   % absorb waves on the top boundary
absorb_bottom=1;% absorb waves on the bottom boundary

%==========================================================================
% path where adjoint sources are located
%==========================================================================

adjoint_source_path='../input/sources/adjoint/';

%==========================================================================
% make wavepropagation movie
%==========================================================================

make_movie='no';                                            % 'yes' or 'no'
movie_file='../output/C_2_heterogeneous_whitened.mp4';     % output file name, should be .mp4