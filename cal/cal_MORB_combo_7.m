% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 10;
cal.nmem   = 15;
cal.nmsy   = 8;
cal.ncmp   = 7;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','Cr$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Cr','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'opx1','opx2','cpx1','cpx2','cpx3','spl1','spl2','for','fay','fsp1','fsp2','fsp3','bio','qtz','phl'};
cal.msyStr = {'opx','cpx','spn','olv','fsp','bio','qtz','phl'};
cal.cmpStr = {'hrz','spx','bas','baf','try','pgr','hab'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 Cr2O3 FeO MgO CaO Na2O K2O H2O  
cal.ioxd = [   1    2    3   4    5   6    7    8   9   10  ]; % oxide indices for viscosity, density functions

%           plg  olv  cpx  qtz
cal.imsy = [ 5    2    1    6];  % mineral system indices for plotting basalt tetrahedron

% oxide composition of mineral end-members
%                SiO2       TiO2      Al2O3      Cr2O3     FeO       MgO        CaO       Na2O      K2O      H2O       
cal.mem_oxd = [  53.53	0.00	7.07	0.00	6.43	30.96	2.01	0.00	0.00	0.00
                   55.41	0.00	3.54	0.00	6.14	31.94	2.97	0.00	0.00	0.00
                   50.77	1.17	0.63	0.00	27.08	5.83	12.88	1.64	0.00	0.00
                   53.24	0.61	6.55	0.00	4.85	23.90	10.06	0.80	0.00	0.00
                   53.83	0.10	2.80	0.00	3.70	19.99	19.58	0.00	0.00	0.00
                    0.00	0.00	41.98	28.04	10.37	19.61	0.00	0.00	0.00	0.00
                    0.00	27.74	0.00	0.00	71.55	0.71	0.00	0.00	0.00	0.00
                   41.07	0.00	0.00	0.00	9.36	49.57	0.00	0.00	0.00	0.00
                   30.87	0.00	0.00	0.00	63.40	5.73	0.00	0.00	0.00	0.00
                   47.76	0.00	33.54	0.00	0.00	0.00	16.53	2.04	0.13	0.00
                   64.14	0.00	22.41	0.00	0.00	0.00	3.52	9.33	0.60	0.00
                   66.11	0.00	19.73	0.00	0.00	0.00	0.85	5.70	7.61	0.00
                   35.95	3.40	11.69	0.00	30.34	5.98	0.00	0.00	9.70	2.94
                  100.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00
                   39.00	2.00	11.00	1.00	4.00	24.00	0.00	0.00	11.00	8.00];           
        


cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100; 

% mineral end-members in mineral systems
% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0  0  0      % orthopyroxene(opx) 
               0  0  1  1  1  0  0  0  0  0  0  0  0  0  0      % clinopyroxene (cpx)
               0  0  0  0  0  1  1  0  0  0  0  0  0  0  0      % spinel (spn) 
               0  0  0  0  0  0  0  1  1  0  0  0  0  0  0      % olivine (olv)
               0  0  0  0  0  0  0  0  0  1  1  1  0  0  0      % feldspar (fsp)
               0  0  0  0  0  0  0  0  0  0  0  0  1  0  0      % biotite (bio)
               0  0  0  0  0  0  0  0  0  0  0  0  0  1  0      % quartz (qtz)
               0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  ];  % phlogopite (phl) 


% mineral end-member composition of melting model components

%Need to populate this table once I have finished MCMC

%              cpx1  cpx2   cpx3    for   fay  spl1  spl2  opx1  opx2  fsp1  fsp2  fsp3  qtz   phl   wat  
cal.cmp_mem = [  0     0      0      0     0     0     0     0     0     0     0     0     0     0    0        % Harzburgite (hrz)
               
                 0     0      0      0     0     0     0     0     0     0     0     0     0     0    0        % Spinel pyroxenite (spx)
                 
                 0     0      0      0     0     0     0     0     0     0     0     0     0     0    0       % Basalt (bas
               
                 0     0      0      0     0     0     0     0     0     0     0     0     0     0    0      % Basalt (baf) another one ...
                 
                 0     0      0      0     0     0     0     0     0     0     0     0     0     0    0      % Trachyte (try)

                 0     0      0      0     0     0     0     0     0     0     0     0     0     0    0     % Plagiogranite (pgr)
               
                 0     0      0      0     0     0     0     0     0     0     0     0     0     0    0 ];    % Hydrous alkaline basalt (hab)
                  
cal.cmp_mem = cal.cmp_mem./sum(cal.cmp_mem,2)*100;


% mineral systems composition of melting model components
cal.cmp_msy = cal.cmp_mem*cal.msy_mem.';

% oxide composition of melting model components
cal.cmp_oxd = cal.cmp_mem*cal.mem_oxd./100;

% oxide composition of mineral systems in melting model components
for i=1:cal.ncmp
    for j=1:cal.nmsy
        cal.cmp_msy_oxd(i,j,:) = cal.cmp_mem(i,cal.msy_mem(j,:)==1)*cal.mem_oxd(cal.msy_mem(j,:)==1,:)./sum(cal.cmp_mem(i,cal.msy_mem(j,:)==1)+1e-32);
    end
end

% set pure component melting points T_m^i at P=0
cal.T0  = [1850 17050 1025 825 800  800750];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [7.5 5.4 2.4 2.0 2.8 2.8 3.0 ];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [7.5 5.4 2.4 2.0 2.8  2.8 3.0];

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [37.0  8.0  6.0  8.0 10.0 10.0  12.0 ];

% set entropy gain of fusion DeltaS [J/K]
cal.Dsx = 300;

% specify geochemical model parameters
cal.ntrc     = 6;                    % number of trace elements
cal.trcStr   = {'K 0.01','K 0.10','K 0.5','K 3.00','K 10.0','K 0.5'};
cal.Ktrc_mem = [0.01;0.10;0.5;3.0;10.0;0.5].*ones(cal.ntrc,cal.nmem);
cal.ElRatio  = [1,4;2,5;3,6];

% specify density parameters
%               for  fay  ant  san  dps  aug  ulv  crm  phl  dol  qtz
cal.rhox0   = [3200,4000,2680,2550,3210,3520,3930,4730,2700,2700,2540]; % mem ref densities [kg/m3]

% specify three-phase coefficient model parameters
%              for  fay  ant  alb   dps  aug  ulv  crm  phl  dol  qtz
cal.etax0   = [1e19,1e19,1e17,1e17,1e20,1e20,1e17,1e17,1e15,1e16,1e19]; % mem ref viscosities [Pas]
cal.Eax     = 300e3;                  % solid viscosity activation energy [J/mol]
cal.AA      =[ 0.5989, 0.1772; ...    % permission slopes
               0.0397, 0.1182 ];      % increases permission slopes away from step function 

cal.BB      =[ 0.6870, 0.3130; ...  % permission step locations
               0.9998, 0.0002;];        % sets midpoint of step functions

cal.CC      =[[0.9826, 0.0174]*9.1697; ... % permission step widths
              [0.1695, 0.8305]*4.2773;];   % factor increases width of step functions

% convergence tolerance
cal.tol     = 1e-9;
cal.alpha   = 0.5;