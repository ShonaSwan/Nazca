% prepare workspace
clear all; close all;

addpath(genpath('../'));
addpath('../../cal')
addpath('../../src')
addpath('../../../unmix')
addpath('../../../unmix/src')
load ocean
Fs = {'FontSize',12};
FS = {'FontSize',15};
FL = {'FontSize',18};
Ms = {'MarkerSize',8};
MS = {'MarkerSize',10};
ML = {'MarkerSize',12};
TX = {'Interpreter','latex'};
TL = {'TickLabelInterpreter','latex'};
LB = {'Location','best'};
LO = {'Location','bestoutside'};
%% *****  prepare for pseudo-component calibration  *********************** （ Crucial！！！）
% A parameter fitting algorithm is used to optimise the result, so that one endmember better represents the differentiation end point, while the other endmember better represents the starting point.

% !!!  Run Section to load end-member calibration and prepare for pseudo-component calibration  !!!

Datamode   =  0; 

if Datamode == 0       % Both trends combined 
    load('MORB_FM_FC_Final');
    cal_MORB_combo_6;  
    
elseif Datamode == 1   % Fractional melting 
   load('MORB_FM_Final');
   cal_MORB_frac_melt;  

elseif Datamode == 2   % Fractional crystalisation  
    load('MORB_FC_Final'); 
    cal_MORB_frac_cryst;  
end 

% !!!  Edit end-member appearances in pseudo-components  !!!
% - number and sequence of end-members must correspond to list in cal file
% - in sequence of appearance, add one new mineral
%   system to next pseudo-component
% - in sequence of appearance, add one new end-member of each mineral
%   system to next pseudo-component
% - phase out mineral systems and their end-members in accordance with
%   their fading or disappearance in PHS_frc


%% Defining the Solidus and Liquidus values for the inital KLB1 composition and the Inital fractional melting composition 

 % KLB-1 mantle: example high-pressure range /% Initial melt: example lower-pressure range
 PslKLB = [3:1:10, 15:5:50]'; PslFC   = ([2,3])';  % kbar 
 % Obtain these temperatures from MAGEMin at each pressure
 TsolKLB = [1139.24,1151.37,1160.84,1167.95,1172.45,1173.73,1170.75,1190.82,1301.46,1388.83,1457.63,1512.81,1557.59,1594.07,1622.67,1644.58 ];
 TliqKLB = [1802.87,1807.59,1812.13,1816.5,1820.7,1824.75,1828.63,1832.38,1849.11,1862.94,1874.35,1883.71,1891.33,1897.43,1902.21,1905.83];
 TsolMelt = [1065.35, 1065.35];
 TliqMelt = [1284.98, 1284.98];
 % TsolMelt = [907.14,909.76,912.47,915.09,917.73,920.37,922.99,925.53,928.12,930.65,933.23];
 % TliqMelt = [1288.33, 1288.45, 1288.57, 1288.67, 1288.79,1288.89,1289.0,1289.11,1289.21,1289.31, 1289.41];

%%
if Datamode == 0 
%                   Opx   |  Cpx   |  Spl  |  Ol  |   Pl  | Bio  |Qrt | Phl
indmem  = logical([ 1 0 0   0 0 0 0   0 0    1 1    0 0 0    0     0    0  ;   % Harzburgite (hrz)
                    1 1 1   1 1 0 0   1 0    1 1    0 0 0    0     0    0  ;   % Spinel pyroxenite (spx)
                    0 0 1   1 1 1 0   1 1    1 1    1 1 0    0     0    0  ;   % Basalt (bas)
                    0 0 0   0 0 1 1   1 1    1 0    0 1 1    0     0    0  ;   % Trachyte (try)
                    0 0 0   0 0 1 1   0 1    0 0    0 1 1    1     1    0  ;   % Plagiogranite (pgr)
                    1 0 0   1 1 0 0   1 1    1 1    1 0 0    0     0    1]);   % Hydrous alkaline basalt (hab)

elseif Datamode == 1 

    %                 Opx |  Cpx   | Spl  |  Ol  |   Phl | Pl
indmem  = logical([ 1 0 0   1 0 0    0     1 1       0   0 0  ;   % Harzburgite (hrz)
                    1 1 1   1 1 0    1     1 1       0   0 0  ;   % Spinel pyroxenite (spx)
                    0 0 1   1 1 1    1     1 1       0   1 1  ;   % Basalt (bas)
                    1 0 0   1 1 0    1     1 1       1   0 0]);   % Hydrous alkaline basalt (hab)
                    
elseif Datamode == 2 
   
    %                 Ol |   Pl   |  Cpx  |  Spl |  Bio  | Qrt 
indmem  = logical([  1 1   1 1 0     1 0    1 1       0      0     ;   % Basalt (bas)
                     1 0   0 1 1     1 1    1 1       0      0     ;   % Trachyte (try)
                     0 0   0 1 1     0 1    0 1       1      1   ]);   % Plagiogranite (pgr)
 
end 

cmp_oxd = 0.5*EMInt + 0.5*EMExt; %*cal.mem_oxd(1:end-1,1:end-1)/100; EM_Int: internal endmembers obtained earlier via PCA and endmember extraction.
%cmp_oxd = [cmp_oxd,zeros(cal.ncmp)];
%cmp_oxd = [cmp_oxd;zeros(1,cal.noxd)];

cmp_mem = zeros(ncmp,cal.nmem);
for ic = 1:ncmp
    imem  = indmem(ic,:);
    A     = cal.mem_oxd(imem,:).';
    b     = squeeze(cmp_oxd(ic,:));
    cmp_mem(ic,imem) = lsqregcmp(A,b,[0.5,5,0,0.002])*100;
end

if Datamode == 0 

    
 %             opx1     opx2      opx3     cpx1      cpx2      cpx3       cpx4      spl1      spl2      for      fay       fsp1      fsp2     fsp3      bio      qtz        phl
 cmp_mem = [22.6000         0         0         0         0         0         0    0.1000         0   76.0000    1.3000         0         0         0         0         0         0
             0.1000   23.2000    0.1000   19.0000    6.7000         0         0    1.8000         0   48.5000    0.6000         0         0         0         0         0         0
                  0         0   16.0400         0         0   25.9300    4.1200         0         0    1.8100    8.7600   38.2700    5.0700         0         0         0         0
                  0         0         0         0         0         0   36.6000    4.6000   12.9000         0         0         0   45.7000    0.2000         0         0         0
                  0         0         0         0         0         0    5.3100         0    0.4000         0         0         0         0   66.0600    5.3100   22.9200         0
             0.1000         0         0   16.2200    9.0100         0         0    3.6000         0   27.1300    3.6000   21.7200         0         0         0         0   18.6200];

elseif Datamode == 1 
            
             % opx1   opx2   opx3    cpx1    cpx2    cpx3    spl1    for    fay      phl     pla1   pla2
cmp_mem = [  22.60	0.00	0.00	0.00	0.00	0.00	0.10	76.00	1.30	0.00	0.00	0.00
              0.10	23.20	0.10	19.00	6.70	0.00	1.80	48.50	0.60	0.00	0.00	0.00
              0.00	0.00	16.04	0.00	0.00	25.93	0.00	1.81	8.76	0.00	38.27	5.07
              0.10	0.00	0.00	16.22	9.01	0.00	3.60	27.13	3.60	18.62	21.72	0.00];

elseif Datamode == 2
          %   for    fay     fsp1    fsp2    fsp3    cpx1   cpx2    spl1     spl2   bio      qtz
cmp_mem = [  1.81	8.76	38.27	5.07	0.00	25.93	4.12	0.00	0.00	0.00	0.00
             0.00	0.00	0.00	45.70	0.20	0.00	36.60	4.60	12.90	0.00	0.00
             0.00	0.00	0.00	0.00	66.06	0.00	5.31	0.00	0.40	5.31	22.92];
                    

end 

%PHS_frc(:,end+1) = 0; 

cmp_mem_init = round(cmp_mem,1);
cmp_mem_init = cmp_mem_init./sum(cmp_mem_init,2)*100;
indmem = logical(cmp_mem_init);
cmp_mem_best = cmp_mem_init;

cmp_oxd_init = cmp_mem_init*cal.mem_oxd/100;
cmp_oxd_best = cmp_oxd_init;


% %Set initial component–endmember conditions for fitting melting point parameters.
if Datamode == 0 
       %  (hrz)  (spx)  (bas)  (try)   (pgr)   (hab)
T0_init = [1850  1600    1120    950     760    700];  T0_best = T0_init; % T0_init = [1850   1663    1180   807  600] 
A_init  = [ 7.5  7.0     2.90    2.00    2.00   4.0];  A_best =  A_init;  %  A_init  = [6.9    5.7   4.5   1.6    3] 
B_init  = [ 7.   7.0     2.50    2.20    1.70   4.6];  B_best =  B_init;  %  B_init  = [3.1 2.8 2.7 2.9 2.8]
r_init  = [ 30    10     1.80    7.7      12    10.0];  r_best =  r_init; %   r_init  = [20.5  11.1   7.9   27.5   7.3]
elseif Datamode == 1 
        %  (hrz)     (spx)    (bas)     (hab)
T0_init = [1850    1600    1120     700];  T0_best = T0_init; 
A_init  = [7.50     7.0    2.90    4.0];  A_best =  A_init;  
B_init  = [7.5      7.0    2.50    4.6];  B_best =  B_init; 
r_init  = [30       10     1.80    10.0];  r_best =  r_init; 

elseif Datamode == 2
       %   (bas)    (try)   (pgr) 
T0_init = [ 1120     995     760];  T0_best = T0_init; % T0_init = [1850   1663    1180   807  600] 
A_init  = [ 2.90    2.00    2.00];  A_best =  A_init;  %  A_init  = [6.9    5.7   4.5   1.6    3] 
B_init  = [ 2.50    2.20    1.70];  B_best =  B_init;  %  B_init  = [3.1 2.8 2.7 2.9 2.8]
r_init  = [ 1.80    7.7       12];  r_best =  r_init; %   r_init  = [20.5  11.1   7.9   27.5   7.3]
end 

% compose initial parameter guess
m0     = [T0_init.';A_init.';B_init.';r_init.';cmp_mem_init(:).*indmem(:);];

% set function to calculate forward model
% m --> dhat
% dhatFunc  = @(model) OxdFromCmpMem(model,MLTp,SOLp,PHS(:,1),cal);
dhatFunc  = @(model) ModelFitP(model,Tmp,Prs,SYS_oxdp,PHS_frc,PslKLB,PslFC,cal,[0.5,6,1.75,1e-4],Datamode); %0.5,5,1.75,1e-3 %

% test fit function for initial guess
[~,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit, PHS_oxdfit,PHS_frcfit,SOL_cmpfit,MLT_cmpfit,SYS_cmpfit,SYS_cmpfit_Anh,TsolfitKLB,TliqfitKLB,TsolfitFC,TliqfitFC,TmKLB,TmFC] = dhatFunc(m0);

%fractional melting 
%evaluate melting points as function of pressure psl = magemin points 

    cal.T0 = T0_init;
    cal.A  = A_init;
    cal.B  = B_init;
    cal.r  = r_init;

if Datamode == 0 
    %fractional melting 
    PP_KLB = linspace(0.001,max(PslKLB),100).';
    var.m = ones(size(PP_KLB))/2;
    var.x = var.m;
    var.c = repmat(SYS_cmpfit_Anh(1,:)./sum(SYS_cmpfit_Anh(1,:)),length(PP_KLB),1);
    var.P = PP_KLB/10;
    var.T = 1000 + PP_KLB*5e-8;
    [~,cal] = leappart(var,cal,'T');
    TmKLB = cal.Tm;

    %fractional crystalstions 
    PP_FC = linspace(min(PslFC),max(PslFC),100).';
    var.m = ones(size(PP_FC))/2;
    var.x = var.m;
    var.c = repmat(SYS_cmpfit_Anh(33,:)./sum(SYS_cmpfit_Anh(33,:)),length(PP_FC),1);
    var.P = PP_FC/10;
    var.T = 1000 + PP_FC*5e-8;
    [~,cal] = leappart(var,cal,'T');
    TmFC = cal.Tm;

elseif Datamode == 1 
    %fractional melting 
    PP_KLB = linspace(0.001,max(PslKLB),100).';
    var.m = ones(size(PP_KLB))/2;
    var.x = var.m;
    var.c = repmat(SYS_cmpfit_Anh(1,:)./sum(SYS_cmpfit_Anh(1,:)),length(PP_KLB),1);
    var.P = PP_KLB/10;
    var.T = 1000 + PP_KLB*5e-8;
    [~,cal] = leappart(var,cal,'T');
    TmKLB = cal.Tm;

elseif Datamode == 2
     %fractional crystalstions 
    PP_FC = linspace(min(PslFC),max(PslFC),100).';
    var.m = ones(size(PP_FC))/2;
    var.x = var.m;
    var.c = repmat(SYS_cmpfit_Anh(1,:)./sum(SYS_cmpfit_Anh(1,:)),length(PP_FC),1);
    var.P = PP_FC/10;
    var.T = 1000 + PP_FC*5e-8;
    [~,cal] = leappart(var,cal,'T');
    TmFC = cal.Tm;

end


% rescale mineral phase fractions
PHS_frc(:,2:end) = PHS_frc(:,2:end)./(100-PHS_frc(:,1)+eps)*100;

% plot basic information for initial fit
level = 3;
run('MCMC/PlotFitFrac.m')

% Fitting algorithm determines the fraction of each mineral endmember in every pseudo component.
%% *****  calibrate pseudo-components and melting point parameters  *******   Monte Carlo algorithm (brun-in & annealing）

% cal_MORB_hi;  % read calibration file
Psl = 1.25;

% uncomment following lines to run MCMC again with previous best fit as initial guess
T0_init = T0_best; A_init = A_best; B_init = B_best; r_init = r_best; cmp_mem_init = cmp_mem_best;
m0      = [T0_init.';A_init.';B_init.';r_init.';cmp_mem_init(:).*indmem(:)];

% *** !!!  set MCMC parameters then Run Section dbstack   to execute MCMC routine  !!!
Niter           = 1e3 ;% 1e4;               % 1e6 number of samples to take; how many samples to take during one Monte Carlo run (how many times you do a random shuffle)
anneal.initstep =   5e-3 ;%1e-3;              % Adjust the step size to achieve a reasonable acceptance ratio of 20–30 per cent. The step size determines how far from the current solution guess the next random perturbation is applied
anneal.levels   = 1;                 % select number of annealing levels
anneal.burnin   = max(1,Niter/1e3);  % Set the length of the initial burn-in sequence. This defines the duration of the burn-in period.
anneal.refine   = max(1,Niter/1e3);  % set length of final refinement sequence

% !!!  set data uncertainties to weight likelihood function  !!!
MLT_scl   = max(0.01,(MLT_oxdp(:)-min(MLT_oxdp(:)))./(max(MLT_oxdp(:))-min(MLT_oxdp(:))));
SOL_scl   = max(0.01,(SOL_oxdp(:)-min(SOL_oxdp(:)))./(max(SOL_oxdp(:))-min(SOL_oxdp(:))));
MEM_scl   = max(0.01,(SOL_mem (:)-min(SOL_mem (:)))./(max(SOL_mem (:))-min(SOL_mem (:))));
PHS_scl   = max(0.01,(PHS_frc (:)-min(PHS_frc (:)))./(max(PHS_frc (:))-min(PHS_frc (:))));
sigma_MLT =  0.1  * MLT_scl.^0.25;       % uncertainty of melt oxide composition
sigma_SOL =  0.1  * SOL_scl.^0.25;       % uncertainty of melt oxide composition
sigma_MEM =  0.1  * MEM_scl.^0.25;       % uncertainty of solid end-member composition
sigma_PHS =  0.1  * PHS_scl.^0.25;       % uncertainty of phase fractions

sigma_TSL = 0.1 * ones(numel([TsolKLB(:); TliqKLB(:); TsolMelt(:); TliqMelt(:)]), 1);
sigma = [sigma_MLT(:); sigma_SOL(:); sigma_MEM(:); sigma_PHS(:); sigma_TSL];
%sigma = [sigma_MLT;sigma_SOL;sigma_MEM;sigma_PHS]; % combine all as in data vector

% load calibration data constraints into data vector
data   = [MLT_oxdp(:);SOL_oxdp(:);SOL_mem(:);PHS_frc(:);TsolKLB(:);TliqKLB(:);TsolMelt(:);TliqMelt(:)];
%data   = [MLT_oxdp(:);SOL_oxdp(:);SOL_mem(:);PHS_frc(:)];

% construct lower and upper parameter bounds
dm    =[1*max(20.0,0.05*T0_init.'); ...
        1*max(0.25,0.4* A_init.'); ...   
        1*max(0.25,0.4* B_init.'); ...
        1*max(0.50,0.2* r_init.'); ... %0.4 
        1*max(5.00,min(20.0,cmp_mem_init(:))).*indmem(:)];
m0_lw  = m0 - dm;
m0_up  = m0 + dm;
mbnds  = [m0_lw(:),m0_up(:)]; % model parameter bounds
mbnds(1*(cal.ncmp)+(1:cal.ncmp)       ,:) = max(0.5,                  mbnds(1*(cal.ncmp)+(1:cal.ncmp)       ,:));
mbnds(2*(cal.ncmp)+(1:cal.ncmp)       ,:) = max(1.0,                  mbnds(2*(cal.ncmp)+(1:cal.ncmp)       ,:));
mbnds(3*(cal.ncmp)+(1:cal.ncmp)       ,:) = max(3.0,                  mbnds(3*(cal.ncmp)+(1:cal.ncmp)       ,:));
mbnds(4*(cal.ncmp)+(1:cal.ncmp*cal.nmem),:) = max(indmem(:)/10,min(99.9,mbnds(4*(cal.ncmp)+(1:cal.ncmp*cal.nmem),:))); 
% The constraint equals 5 times the number of components, minus 1 (to handle cases like 1, 2, 3, 4, 5), and then multiplied by a parameter for each component except water. The number of components times the number of endmembers gives the edges in the table. If an initial value (from the earlier initialization) falls outside the distribution boundary (=0), the MCMC can fail or get stuck in this exact way.

mbnds(m0==100 ,:) = 100;
% mbnds(m0==7.0 ,1) = 5.0;


mbnds(m0==T0_init(1),:) = T0_init(1);
mbnds(m0==T0_init(end),:) = T0_init(end);
anneal.initstep = anneal.initstep * diff(mbnds,1,2);  % resize step according to bounded bracket

% set parameter names according to info from calibration file
mNames = cell(cal.ncmp*cal.nmem+4*(cal.ncmp),1);
k = 1;
for j=1:4
    for i=1:cal.ncmp
        if j==1
            mNames{k} = ['T0:',cal.cmpStr{i}];
        elseif j==2
            mNames{k} = ['A:',cal.cmpStr{i}];
        elseif j==3
            mNames{k} = ['B:',cal.cmpStr{i}];
        elseif j==4
            mNames{k} = ['r:',cal.cmpStr{i}];
        end
        k = k+1;
    end
end
for j=1:cal.nmem
    for i=1:cal.ncmp
        mNames{k} = [cal.memStr{j},':',cal.cmpStr{i}];
        k = k+1;
    end
end

% set function to apply further constraints to a proposed set of parameter values
% m --> m
ConstrFunc = @(model) ConstrFuncs('SumConstr', model, cal.ncmp, cal.nmem, 100);

% set function to calculate prior probability given a set of model param values
% m --> prior prob
PriorFunc = @(model) ProbFuncs('PriorFunc', model, mbnds, 'uniform');

% set function to calculate likelihood of forward model
% dhat --> likelihood 
LikeFunc  = @(dhat,model) ProbFuncs('LikeFuncSimplex',dhat,data,sigma,0.1,0.1,max(Psl)*3,model,cal);

bestfit = m0;  % initialise bestfit from initial conditions


%*****  RUN MCMC PARAMETER FITTING ROUTINE  *******************************
% This line passes the burn-in and annealing information to the MCMC routine. The MCMC reads anneal.burnin and treats the first Niter/5 iterations as the burn-in period (usually with a more relaxed acceptance rate and larger step size, as implemented internally in MCMC).
tic;
[models,prob,accept,bestfit] = mcmc(dhatFunc,PriorFunc,LikeFunc,ConstrFunc,4*(cal.ncmp)+(1:cal.ncmp*cal.nmem),m0,mbnds,anneal,Niter);
RunTime(1) = toc;
%**************************************************************************


% uncomment following line to plot likelihood histograms for fitted parameters (slow!)
% plotmcmc(models, prob, [], mbnds, anneal, mNames); 

% extract best fit parameters
T0_best       = bestfit(               (1:cal.ncmp)).';
A_best        = bestfit(1*(cal.ncmp)+(1:cal.ncmp)).';
B_best        = bestfit(2*(cal.ncmp)+(1:cal.ncmp)).'; %A_best;
r_best        = bestfit(3*(cal.ncmp)+(1:cal.ncmp)).';  % a set of fitting parameters which have to do with the curvature of the phase loops that this model is generating.
cmp_mem_best  = reshape(bestfit(4*(cal.ncmp)+(1:cal.ncmp*cal.nmem)),cal.ncmp,cal.nmem);
cmp_oxd_best  = cmp_mem_best*cal.mem_oxd/100;

% evaluate forward model for best fit parameters
[dhat,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmpfit,MLT_cmpfit,SYS_cmpfit,TsolfitKLB,TliqfitKLB,TsolfitFC,TliqfitFC,TmKLB,TmFC] = dhatFunc(bestfit);
%[dhat,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmpfit,MLT_cmpfit,SYS_cmpfit] = dhatFunc(bestfit);
[Lbest,Vsimplex] = LikeFunc(dhat,bestfit);

% evaluate melting points as function of pressure
% if isfield(cal,'Tsol'); cal = rmfield(cal,{'Tsol' 'Tliq'}); end
% PP         = linspace(0.001,max(Psl)*2,50).';  
% var.m      = ones(size(PP))/2; var.x = var.m; var.f = 0*var.m;
% cal.T0     = T0_best;
% cal.A      = A_best;
% cal.B      = B_best;
% cal.r      = r_best;
% cal.dTH2O  = dT_best;
% var.c      =repmat(SYS_cmpfit(1,:).*[ones(1,cal.ncmp-1),0]./sum(SYS_cmpfit(1,1:end-1),2),length(PP),1);   % component fractions [wt] 
% var.P      = PP/10;
% var.T      = 1000+PP*5e-8;
% [~,cal]    = leappart(var,cal,'T');
% Tm         = cal.Tm;

% uncomment following line to retrieve probability distributions (slow!)
% Nbins = min(500,Niter/20);
% [ppd_mcmc.m, ppd_mcmc.prob] = CalcPDF(mbnds, models(anneal.burnin:end,:), Nbins);
% 

%% *****  visualise best fit calibration  *********************************

%!!!  adjust desired level of detail to plot then run section  !!!
%     level = 1     only simple line plots
%     level = 2     add Harker diagrams and T-X pseudo-sections
%     level = 3     add mineral systems and display best fit parameters

level = 3;
run('MCMC/PlotFit.m')

%% save and display calibration
save('MtAp_noamp10_calibration_2');