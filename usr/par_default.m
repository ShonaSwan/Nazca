% set run parameters
runID    =  'default';           % run identifier
srcdir   =  '../src';            % output directory
outdir   =  '../out';            % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence
colourmap = 'ocean';             % choose colourmap ('ocean','acton','devon','lajolla','lipari','lapaz','glasgow')

% set model domain parameters
D         =  200e3;               % System depth [m]
L         =  1*D;                 % System width [m]
N         =  100;                 % number of grid points in z-direction (incl. 2 ghosts)
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
Nt       =  5e5;                 % number of time steps to take
hr       =  3600;                % conversion seconds to hours
yr       =  24*365.25*hr;        % conversion seconds to years
tend     =  1e9*yr;              % end time for simulation [s]
dt       =  1e2*yr;              % initial time step [s]
dtmax    =  1e6*yr;              % maximum time step [s]
minage   =  20e6*yr;             % age of system before simuation 

% set up melt fraction variabels 
minit     =  0.01;                % maximum initial melt fraction (Initial reduction of melt)
mumin     =  1e-5;                % Setting lower limit for melt fraction in 
mumax     =  0.25;                % Setting upper limit for melt fraction in 

% set up mid ocean ridge spreading 
bnd_sprc  =  6e3;                 % Top boundary horizontal coordinate (centre) of spreading rate 'S' function [km]  
bnd_sprw  =  5e3;                 % Width of top boundary spreading rate 'S' function [km] 
sprate   =  0.04/yr;             % Half spreading rate [m/s] (modeling half the ridge)

% set initial thermo-chemical state (Mantle)
init_mode = 'plume';              % 'plume' or 'MOR'
bndmode   =  1;                   % boundary assimilation mode (0 = MOR; 1 = Plume)
seed      =  24;                  % random perturbation seed
smth      =  10;                  % regularisation of initial random perturbation
zlay      =  2.0;                 % layer thickness (relative to domain depth D)
dlay      =  0.0;                 % random perturbation to layer thickness (relative to grid spacing h)
wlay_T    =  0*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
wlay_c    =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
T0        =  5;                   % temperature top layer [deg C]
T1        =  1350;                % temperature base layer [deg C]
dTr       =  0;                   % amplitude of random noise [deg C]
dTg       =  0;                   % amplitude of centred gaussian [deg C]
c0        =  [0.73 0.16 0.10 0.01 0];  % components (maj comp, H2O) top layer [wt] (will be normalised to unit sum!)
c1        =  c0;                  % components (maj comp, H2O) bot layer [wt] (will be normalised to unit sum!)
dcr       =  [1,-1,0,0,0]*0e-3;     % amplitude of random noise [wt SiO2]
dcg       =  [0,0,0,0,0,0,0];     % amplitude of centred gaussian [wt SiO2]

% set model trace and isotope geochemistry parameters (must match # trace elements and isotope ratios in calibration!)
trc0     =  [1,1,1,1,1,1];       % trace elements top layer [wt ppm]
trc1     =  [1,1,1,1,1,1];       % trace elements base layer [wt ppm]
dr_trc   =  [0,0,0,0,0,0];       % trace elements random noise [wt ppm]
dg_trc   =  [0,0,0,0,0,0];       % trace elements centred gaussian [wt ppm]

% set initial thermo-chemical state (Crust)
crust_sw  =  1;                     % 0 = no crust, 1 = crust 
Hcmin     =  6e3;                   % Minimum crustal thickness 
c_crust   =  [0.01 0.13 0.80 0.06 0];    % components (maj comp, H2O) Crustal layer
trc_crust =  [0.1,0.1,0.5,10,10,2]; % trace elements crust layer [wt ppm]

% set initial thermo-chemical state (Plume)
dT_plume  = 150;                                % Temperature difference between the plume and the mantle 
pl_width  = 50e3;                               % Width of the plume [m]
pl_local  = L/2;                                % Location of the mantle plume along the bottom boundary [m]
c_plume   = [0.70 0.18 0.11 0.01 0];            % components of plume (maj comp, H2O) [wt] (will be normalised to unit sum!)
trc_plume = [10.0, 10.0, 2.0, 0.1, 0.1, 2.0];   % trace elements system plume [wt ppm]

% Plastic Deformation 
tyield    =  5e7;                 % yield stress for shear failure [Pa]
pyield    =  1e7;                 % yield pressure for tensile failure [Pa]
etaymin   =  1e20;                % minimum yield viscosity

% Melt Extraction Algorythm 
erupt_ratio = 0.5;                % 1 = all eruption (surface), 0 = all 
tau_e     =  1e4*yr;              % extraction/eruption time (set to 0 to tie to dt)
mthr      =  0.10;                % threshold melt fraction for extraction/eruption

% set thermo-chemical boundary parameters
% fractxtl =  0;                   % fractional crystallisation mode for 0-D (Nz=Nx=1)
% fractmlt =  0;                   % fractional melting mode for 0-D (Nz=Nx=1)
% fractres =  0.25;                % residual fraction for fractionation mode
% dPdT     =  0e5;                 % decompression rate for 0D models

Ptop     =  4.0e7;                 % top pressure [Pa]

bnd_w    =  h/2;                   % boundary layer width [m]
bnd_h    =  [0,0,0];               % internal wall rock layer thickness [m]
tau_T    =  1e4*yr;                % wall cooling/assimilation time [s]
tau_a    =  24*hr;                 % wall cooling/assimilation tie [s]
Twall    =  [T0,nan,nan,nan];      % [top,bot,left,right] wall rock temperature [degC] (nan = insulating)
cwall    =  nan(3,7,7);            % [top,bot,left,right] wall rock major component [wt SiO2] (nan = no assimilation)
trcwall  =  nan(3,6,6);            % [top,bot,left,right] wall rock trace elements [wt ppm] (nan = no assimilation)

% set thermo-chemical material parameters
calID    =  'MORB_lo';           % phase diagram calibration
aTm      =  5e-5;                % melt  thermal expansivity [1/K]
aTx      =  1e-5;                % xtal  thermal expansivity [1/K]
kTm      =  1;                   % melt  thermal conductivity [W/m/K]
kTx      =  5;                   % xtal  thermal conductivity [W/m/K]
cPm      =  1300;                % melt  heat capacity [J/kg/K]
cPx      =  1000;                % xtal  heat capacity [J/kg/K]
tau_r    =  1e3*yr;              % reaction time scale (set to zero for quasi-equilibrium mode)

% set model buoyancy and pressure parameters
bPx      =  1e-11;               % solid compressibility [1/Pa]
bPm      =  1e-11;               % melt  compressibility [1/Pa]
dx0      =  1e-2;                % crystal size [m]
g0       =  10.;                 % gravity [m/s2]

% set numerical model parameters
meansw   =  0;                   % 0 = Geometric mean 1 = Arithmetic mean
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
maxcmp   =  0.01;                % maximum change in phase fraction due to compaction
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-8;                % outer its absolute tolerance
maxit    =  15;                  % maximum outer its
alpha    =  0.40;                % iterative step size parameter
beta     =  0.00;                % iterative damping parameter
gamma    =  0.20;                % iterative lagging of viscosities
lambda1  =  0e-7;                % pressure regularisation parameter
lambda2  =  0e-7;                % pressure regularisation parameter
etacntr  =  1e5;                 % maximum viscosity contrast
Delta_cnv=  h/2;                 % correlation length for eddy diffusivity (multiple of h, 0.5-1)
Delta_sgr=  dx0*10;              % correlation length for phase fluctuation diffusivity (multiple of dx0, df0, 10-20)
Prt      =  3;                   % turbulent Prandtl number (ratio of momentum to heat diffusivity)
Sct      =  3;                   % turbulent Schmidt number (ratio of momentum to mass diffusivity)
etamin   =  1e18;                % minimum viscosity
kmin     =  1e-9;                % minimum diffusivity
kmax     =  1e+9;                % maximum diffusivity
Pcouple  =  0;                   % coupling phase equilibria and material properties to dynamic pressure
Rcouple  =  0;                   % switch on for full reactive coupling

% set various options
calibrt  =  0;                   % not in calibrate mode
bnchm    =  0;                   % not a benchmark run

