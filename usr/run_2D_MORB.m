 % prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID     =  '2D_MORB_N100';      % run identifier
restart   =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop       =  10;                  % output frame plotted/saved every 'nop' time steps
plot_op   =  1;                   % switch on to live plot results
save_op   =  1;                   % switch on to save output to file
plot_cv   =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D         =  200e3;               % chamber depth [m]
N         =  100;                 % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  1.0*D;               % chamber width (equal to h for 1-D mode) [m]
sprate    =  0.03/yr;             % Half spreading rate [m/s] (modeling half the ridge)
Hcmin     =  0e3;  %6e3               % Minimum crustal thickness 
bnd_sprc  =  6e3;                 % Top boundary horizontal coordinate (centre) of spreading rate 'S' function [km]  
bnd_sprw  =  5e3;                 % Width of top boundary spreading rate 'S' function [km] 

% set model timing parameters
Nt        =  5e5;                 % number of time steps to take
tend      =  1e9*yr;              % end time for simulation [s]
dt        =  1e2*yr;              % initial time step [s]
mumin     =  1e-5;                % Setting lower limit for melt fraction in coeff. 
mumax     =  0.2;                 % Setting upper limit for melt fraction in coeff.

% set initial thermo-chemical state
init_mode =  'MOR';
minage    =  7e5*yr;
T0        =  5;                   % temperature top  layer [deg C]
T1        =  1350;                % temperature base layer [deg C]
wlay_c    =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
c0        =  [0.82 0.17 0.01 0];  % components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c_crust   =  [0.01 0.90 0.09 0];  % components (maj comp, H2O) Crustal layer
c1        =  c0;                  % components (maj comp, H2O) base layer [wt] (will be normalised to unit sum!)
dcr       =  [1,-1,0,0]*0e-3;     % Random perturbation of the composition field
dr_trc    =  [0,0,0,0,0,0];       % trace elements random noise

% set model trace and isotope geochemistry parameters (must match # trace elements and isotope ratios in calibration!)
trc0      =  [1,1,1,1,1,1];       % trace elements system layer [wt ppm]
trc_crust =  [0.1,0.1,0.5,10,10,2];       % trace elements crust layer [wt ppm]

% set thermo-chemical boundary parameters
periodic  =  0;
bndmode   =  5;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = Mid-Ocean Ridge setup)
bnd_w     =  h;                   % boundary layer width [m]
tau_T     =  1e4*yr;              % wall cooling/assimilation time [s]
Twall     =  [T0,nan,nan,nan];    % [top,bot,left,right] wall rock temperature [degC] (nan = insulating)
cwall     =  nan(3,7,7);
Ptop      =  4.0e7;               % top pressure [Pa]

% set thermo-chemical material parameters
calID     =  'MORB_lo';           % phase diagram calibration
tau_r     =  0;                   % phase change reaction time (set to 0 to tie to dt)
tau_e     =  0;                   % extraction/eruption time (set to 0 to tie to dt)
mthr      =  0.20;                % threshold melt fraction for extraction/eruption
minit     =  0.01;                % maximum initial melt fraction

% physical parameters
dx0       =  5e-3;                % matrix grain size
aTm       =  5e-5;                % melt  thermal expansivity [1/K]
aTx       =  1e-5;                % xtal  thermal expansivity [1/K]
kTm       =  1;                   % melt  thermal conductivity [W/m/K]
kTx       =  5;                   % xtal  thermal conductivity [W/m/K]
cPm       =  1300;                % melt  heat capacity [J/kg/K]
cPx       =  1000;                % xtal  heat capacity [J/kg/K]
tyield    =  1e8;                 % yield stress for shear failure [Pa]
pyield    =  1e8;                 % yield pressure for tensile failure [Pa]
etaymin   =  1e20;                % minimum yield viscosity

% set numerical model parameters
TINT      =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN      =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL       =  0.75;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol      =  1e-4;                % outer its relative tolerance
atol      =  1e-8;                % outer its absolute tolerance
maxit     =  15;                  % maximum outer its
alpha     =  0.75;                 % iterative step size
gamma     =  0.25;                 % relaxing parameter for viscosity update
etacntr   =  1e5;                 % maximum viscosity contrast
etamin    =  1e18;                % minimum viscosity
Rcouple   =  0;                   % switch on for full reactive coupling
Pcouple   =  0;                   % switch on for full pressure coupling

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

