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
plot_cv   =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D         =  50e3;               % chamber depth [m]
N         =  100;                 % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  D/5;                   % chamber width (equal to h for 1-D mode) [m]
sprate    =  0.0/yr;              % Half spreading rate [m/s] (modeling half the ridge)
bnd_sprc  =  5e3;                 % Top boundary horizontal coordinate (centre) of spreading rate 'S' function [km]  
bnd_sprw  =  2e3;                 % Width of top boundary spreading rate 'S' function [km] 

% set model timing parameters
Nt        =  5e5;                   % number of time steps to take
tend      =  1e9*yr;                % end time for simulation [s]
dt        =  1e3*yr;                % initial time step [s]

% set initial thermo-chemical state
init_mode =  'linear';
m0        =  0.02;
dmr       =  0.00;
dmg       =  0.00;
T0        =  1350;                   % temperature top  layer [deg C]
T1        =  1350;                % temperature base layer [deg C]
c0        =  [0.85 0.15 0];       % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1        =  c0;                  % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
dcr       =  [1,-1,0]*0e-4;       % Random perturbation of the composition field
dr_trc    =  [0,0,1,0,0,-1];      % trace elements random noise
reactive  =  0;                   % 1 for reactive flow, 0 for non-reactive flow (melt model on off switch)

% set thermo-chemical boundary parameters
periodic  =  0;
bndmode   =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = Mid-Ocean Ridge setup)
bnd_w     =  h;                   % boundary layer width [m]
tau_T     =  1e5*yr;              % wall cooling/assimilation time [s]
Twall     =  [nan,nan,nan,nan];    % [top,bot,left,right] wall rock temperature [degC] (nan = insulating)
cwall     =  nan(3,7,7);
Ptop      =  4.0e7;               % top pressure [Pa]
fin       =  0;
fout      =  1;

% set thermo-chemical material parameters
calID     =  'MORB_lo';              % phase diagram calibration

aTm = 0;
aTx = 0;
bPm = 0;
bPx = 0;

% set numerical model parameters
TINT      =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN      =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL       =  1.0;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol      =  1e-4;                % outer its relative tolerance
atol      =  1e-7;                % outer its absolute tolerance
maxit     =  15;                  % maximum outer its
gamma     =  0;
dtmax     =  1e6*yr;
etacntr   =  1e6;
etamin    =  1e17;

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

