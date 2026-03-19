% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 10;
cal.nmem   = 18;
cal.nmsy   = 6;
cal.ncmp   = 5;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','Cr$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Cr','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'for','fay','ant','alb','san','dps','aug','ulv','mgt','ilm','c','b','a','d','gar','garm','qtz','wat'};
cal.msyStr = {'gar','cpx','olv','opx','spl','pla'};
cal.cmpStr = {'dun','tro','gbr','oth','vol'};


for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 Cr2O3 FeO MgO CaO Na2O  H2O
cal.ioxd = [   1    2    3    4    5   6    7   8   9  ]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                  SiO2     TiO2      Al2O3     Cr2O3     FeO        MgO       CaO      Na2O     K2O        H2O
cal.mem_oxd = [ 42.3700    0.6300   21.3900    2.1000    6.3500   21.0900    6.0700         0         0         0
                42.2400    0.6700   20.8100    2.7100    6.5200   20.9900    6.0600         0         0         0
                42.1500    1.0300   20.1800    2.1400    7.6500   21.0400    5.8100         0         0         0
                52.4300         0    7.8000         0    5.2100   22.1800   11.4200    0.9600         0         0
                57.3800         0    2.0500         0    4.7500   26.6500    7.7300    1.4400         0         0
                53.1800    0.0300    2.3600    0.0500    2.9700   19.0500   22.3600         0         0         0
                52.5100    0.8700    2.6900    0.4100    5.9000   14.9000   22.1800    0.5300    0.0100         0
                41.1500         0         0         0    8.6900   50.1600         0         0         0         0
                38.5000         0         0         0   21.7000   39.1200    0.6800         0         0         0
                53.0600         0    7.3700    0.5500    6.2000   30.7300    2.0900         0         0         0
                54.2100         0    4.1900    1.4600    5.8200   31.4400    2.8800         0         0         0
                57.7900         0    1.2700    0.1400    5.2100   34.0800    1.5100         0         0         0
                      0         0   51.4400   18.0100    9.1200   21.4300         0         0         0         0
                      0         0   26.2500   44.6900   11.9600   17.1000         0         0         0         0
                      0    0.1800   33.7900   34.3300   22.9700    8.7300         0         0         0         0
                46.6800         0   34.3000         0         0         0   17.4000    1.6200         0         0
                48.3300         0   33.1900         0         0         0   16.1100    2.3700         0         0
                      0         0         0         0         0         0         0         0         0  100.0000]; % water (wat)

cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100; 

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0 0  0 0 0 0  0   % garnet (gar)
               0  0  1  1  1  0  0  0  0  0  0  0 0  0 0 0 0  0   % clinopyroxene (cpx)
               0  0  0  0  0  1  1  0  0  0  0  0 0  0 0 0 0  0   % olivine (olv)
               0  0  0  0  0  0  0  1  1  1  0  0 0  0 0 0 0  0   % orthopyroxen (opx)
               0  0  0  0  0  0  0  0  0  0  1  0 0  0 0 0 0  0   % spinel (spl)
               0  0  0  0  0  0  0  0  0  0  1  0 0  0 0 0 0  0]; % plagioclase (pla)
            


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
