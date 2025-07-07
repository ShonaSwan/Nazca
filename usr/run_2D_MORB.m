 % prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID     =  '2D_MORB';           % run identifier
restart   =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop       =  1;                  % output frame plotted/saved every 'nop' time steps
plot_op   =  1;                   % switch on to live plot results
save_op   =  1;                   % switch on to save output to file
plot_cv   =  1;                   % switch on to live plot iterative convergence

% set model domain parameters
D         =  200e3;               % chamber depth [m]
N         =  100;                 % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  1*D;                 % chamber width (equal to h for 1-D mode) [m]
sprate    =  0.04/yr;             % Half spreading rate [m/s] (modeling half the ridge)
bnd_sprc  =  6e3;                 % Top boundary horizontal coordinate (centre) of spreading rate 'S' function [km]  
bnd_sprw  =  3e3;                 % Width of top boundary spreading rate 'S' function [km] 

% set model timing parameters
Nt        =  5e5;                   % number of time steps to take
tend      =  1e9*yr;                % end time for simulation [s]
dt        =  1e2*yr;                % initial time step [s]
mulim     =  1e-6;                   % Setting a limint for melt fraction

% set initial thermo-chemical state
init_mode =  'MOR';
T0        =  5;                   % temperature top  layer [deg C]
T1        =  1350;                % temperature base layer [deg C]
c0        =  [0.85 0.15 0];       % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1        =  c0;                  % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
dcr       =  [1,-1,0]*1e-4;       %Random perturbation of the composition field
dr_trc    =  [0,0,0,0,0,0];      % trace elements random noise
reactive  =  1;                   % 1 for reactive flow, 0 for non-reactive flow (melt model on off switch)    

% set thermo-chemical boundary parameters
periodic  =  0;
bndmode   =  5;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = Mid-Ocean Ridge setup)
bnd_w     =  h;                   % boundary layer width [m]
tau_T     =  1e4*yr;              % wall cooling/assimilation time [s]
Twall     =  [T0,nan,nan,nan];    % [top,bot,left,right] wall rock temperature [degC] (nan = insulating)
cwall     =  nan(3,7,7);
Ptop      =  4.0e7;               % top pressure [Pa]
fin       =  0;
fout      =  0;

% set thermo-chemical material parameters
calID     =  'MORB_lo';              % phase diagram calibration

% grain size
dx0 = 3e-3;
dm0 = 3e-3;
aTm = 5e-5;
aTx = 1e-5;
bPm = 3e-11;
bPx = 1e-11;

% set numerical model parameters
TINT      =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN      =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL       =  0.5;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol      =  1e-4;                % outer its relative tolerance
atol      =  1e-7;                % outer its absolute tolerance
maxit     =  15;                  % maximum outer its
gamma     =  0;
dtmax     =  1e6*yr;
etacntr   =  1e6;
etamin    =  1e17;
alpha     =  0.75;

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

