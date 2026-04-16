% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 10;
cal.nmem   = 19;
cal.nmsy   = 6;
cal.ncmp   = 5;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','Cr$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Cr','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'for','fay','ant','alb','san','dps','aug','ulv','mgt','ilm','c','b','a','d','e','gar','garm','qtz','wat'};
cal.msyStr = {'opx','cpx','olv','gar','spl','pla'};
cal.cmpStr = {'dun','tro','gbr','oth','vol'};


for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 Cr2O3 FeO MgO CaO Na2O  H2O
cal.ioxd = [   1    2    3    4    5   6    7   8   9  ]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                  SiO2     TiO2      Al2O3     Cr2O3     FeO        MgO       CaO      Na2O     K2O        H2O
cal.mem_oxd = [ 52.4000         0    8.0700    0.7000    6.5900   30.2400    2.0000         0         0         0
                55.1500         0    2.8000    1.2000    5.9600   31.7200    3.1700         0         0         0
                58.0000         0    0.9600    0.0300    5.2600   34.3400    1.4100         0         0         0
                43.4600         0   18.4500         0    6.3000   16.6100   15.1800         0         0         0
                53.6300    0.2100    3.0000         0    1.3800   21.0400   20.7400         0         0         0
                58.8200         0         0         0    3.1900   25.6600   10.7900    1.5400         0         0
                48.2800    1.1000    1.0600         0   30.7200         0   17.1100    1.7300         0         0
                41.1900         0         0         0    8.4800   50.3300         0         0         0         0
                30.8900         0         0         0   63.4400    5.6700         0         0         0         0
                42.2100    0.6400   21.0700    2.5900    6.5600   20.6100    6.3200         0         0         0
                42.1200    1.0200   20.5200    2.0600    7.4600   20.6800    6.1400         0         0         0
                42.0500    1.0400   20.3700    1.9600    8.0200   20.4400    6.1200         0         0         0
                      0         0   55.2100   19.2800         0   25.5100         0         0         0         0
                      0         0    8.4500   49.4500   42.1000         0         0         0         0         0
                      0   35.9500         0         0   63.0700    0.9800         0         0         0         0
                43.0800         0   36.7000         0         0         0   20.2200         0         0         0
                61.0100         0   24.6400         0         0         0    6.0900    8.2600         0         0
                63.7600         0   22.4200         0         0         0    3.6200    8.5400    1.6600         0
                      0         0         0         0         0         0         0         0         0  100.0000]; % water (wat)

cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100; 

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0 0  0 0 0 0  0 0  % orthopyroxen (opx)
               0  0  1  1  1  0  0  0  0  0  0  0 0  0 0 0 0  0 0 % clinopyroxene (cpx)
               0  0  0  0  0  1  1  0  0  0  0  0 0  0 0 0 0  0 0 % olivine (olv)
               0  0  0  0  0  0  0  1  1  1  0  0 0  0 0 0 0  0 0  % garnet (gar)
               0  0  0  0  0  0  0  0  0  0  1  0 0  0 0 0 0  0 0 % spinel (spl)
               0  0  0  0  0  0  0  0  0  0  1  0 0  0 0 0 0  0 0];   % plagioclase (pla))
              

            
cal.cmp_mem = zeros(cal.ncmp,cal.nmem);
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
cal.T0  = [1890  1180  1161  1081  986  820];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [6.7  5.4  5.1  2.8  2.3  1.4];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [6.5  5.1  4.3  2.7  1.7  1.3];

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [31.0  3.0  3.0  6.8  9.3  6.0];

% set entropy gain of fusion DeltaS [J/K]
cal.Dsx  = 350;
cal.Dsf  = 450;

% specify melting point dependence on H2O
cal.dTH2O = [889  1418  1455  1556  1697  2049];  % solidus shift from water content prefactor [K/wt^pH2O]
cal.pH2O  = 0.75;                                  % solidus shift from water content exponent

% primary and evolved end-member compositions used in calibration
cal.c0     = [0.044  0.248  0.269  0.317  0.100  0.022  0.005];
cal.c1     = [0.001  0.001  0.001  0.001  0.299  0.697  0.024];

cal.c0_oxd = [50.12  1.01  15.09  9.05  10.57  11.38  2.68  0.10  0.50];
cal.c1_oxd = [76.34  0.16  11.84  2.80   0.71   1.61  4.42  2.12  2.40];

% specify geochemical model parameters
cal.ntrc     = 6;                    % number of trace elements
cal.trcStr   = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%               for  fay  ant  alb  san  dps  aug  ulv  mgt  ilm  hyp  fsl  qtz  wat
cal.rhox0   = [3200,4000,2680,2600,2550,3200,3470,3880,4650,4720,3410,3660,2540,1000]; % mem ref densities [kg/m3]
cal.rhof0   = 1000;                 % fluid ref density [kg/m3]

% specify three-phase coefficient model parameters
%              for  fay  ant  alb  san  dps  aug  ulv  mgt  ilm  hyp  fsl  qtz  wat
cal.etax0   = [1e19,1e19,1e17,1e17,1e17,1e19,1e19,1e17,1e17,1e17,1e19,1e19,1e19,1e0]; % mem ref viscosities [Pas]
cal.etaf0   = 0.1;                    % fluid viscosity constant [Pas]
cal.Eax     = 300e3;                  % solid viscosity activation energy [J/mol]
cal.AA      =[ 0.65, 0.25, 0.35; ...  % permission slopes
               0.20, 0.20, 0.20; ...  % generally numbers between 0 and 1
               0.20, 0.20, 0.20; ];   % increases permission slopes away from step function 

cal.BB      =[ 0.55, 0.18, 0.27; ...  % permission step locations
               0.64,0.012,0.348; ...  % each row sums to 1
               0.80, 0.12, 0.08; ];   % sets midpoint of step functions

cal.CC      =[[0.30, 0.30, 0.40]*0.7; ... % permission step widths
              [0.52, 0.40, 0.08]*1.1; ... % square brackets sum to 1, sets angle of step functions
              [0.15, 0.25, 0.60]*0.7; ];  % factor increases width of step functions

% convergence tolerance
cal.tol     = 1e-9;
cal.alpha   = 0.5;
