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
load('MORB_Comb_6');

cal_MORB_combo_6;  % read calibration file

% !!!  Edit end-member appearances in pseudo-components  !!!
% - number and sequence of end-members must correspond to list in cal file
% - in sequence of appearance, add one new mineral
%   system to next pseudo-component
% - in sequence of appearance, add one new end-member of each mineral
%   system to next pseudo-component
% - phase out mineral systems and their end-members in accordance with
%   their fading or disappearance in PHS_frc
% 
%                     Cpx   |  Ol |  Spl  |  Opx  |  Pl  |  Qrt | Phl
indmem  = logical([ 1 0 0    1 1    0 0      1 0    0 0 0    0     0  0  ;   % dunite (dun)
                    1 1 0    1 1    1 0      1 1    0 0 0    0     0  0  ;   % spinel pyroxenite (spx) 
                    0 1 1    1 1    1 1      0 0    1 1 0    0     0  0  ;   % basalt (bas)
                    0 1 1    0 1    0 1      0 0    1 1 1    0     0  0  ;   % basalt (bas) evolved
                    0 0 1    0 0    0 1      0 0    0 1 1    1     0  1  ;   % plagiogranite (pgr)
                    1 0 1    1 1    1 0      1 0    1 0 0    0     1  1]);   % hydrous alkali basalt (hab)
               
% %                     Cpx   |  Ol |  Spl  |  Opx  |  Pl  |  Qrt |   Phl
% indmem  = logical([ 1 0 0    1 1    1 0     1 0     0 0 0    0     0  0  ;   % dunite (dun)
%                     1 1 1    1 1    1 0     1 1     0 0 0    0     0  0  ;   % spinel pyroxenite (spx) 
%                     0 1 1    0 1    0 1     0 0     1 1 0    0     0  0  ;   % basalt (bas)
%                     0 0 1    0 0    0 1     0 0     1 1 1    1     0  1  ;   % plagiogranite (pgr)
%                     1 0 0    1 1    1 1     1 0     1 0     0     1  1]);   % hydrous alkali basalt (hab)


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

% % getting the fractional melting to work maybe need garnet pyroxenite ?

%            %cpx1      cpx2      cpx3        for       fay        spl1     spl2    opx1     opx2      fsp1      fsp2      fsp3       qtz       phl       wat 
% cmp_mem = [       0         0         0   88.0000   12.0000         0         0         0         0         0         0         0         0         0         0
%             20.6000    0.1000    0.2000   67.1000    5.0000    1.7000         0    9.1000    3.2000         0         0         0         0         0         0
%                   0   40.1000   15.8000    3.9000    2.5000    0.1000    0.1000         0         0   31.4000    0.1000         0         0         0         0
%                   0         0   32.4000         0         0         0    6.1000         0         0    0.1000   30.9000   40.8000    25.7000        0     5.1000
%             42.1000         0         0    0.2000    5.9000    0.6000         0    6.5000         0   37.8000         0         0         0    10.6000    2.3000];

%            %cpx1      cpx2      cpx3        for       fay        spl1     spl2    opx1     opx2      fsp1      fsp2      fsp3       qtz       phl       wat 
% cmp_mem = [       0         0         0   92.0000   8.0000         0        0         0         0         0         0         0         0         0         0
%                25.0         0         0   15.1000   0              2        1        35        17         0         0         0         0         0         0
%                   0   35.1000        10         5          3      0.5     0.5          0         0   35.4000    10.1000         0         0         0         1
%                   0         0         2         0         0         0    6.1000         0        0    0         35.1000   33.000   28.8000            0     2
%             35.1000         0         0    10.2000    5.9000    1.6000         0    0         0   39.8000         0         0         0    5    5];


% %                     cpx1       cpx2      cpx3    cpx4       for        fay      spl1       spl2     spl3     opx1       opx2      opx3    fsp1      fsp2      fsp3       qtz       phl       wat      
% cmp_mem = [         10.4000         0         0         0   77.7000   15.5000         0         0         0         0    8.5000         0         0         0         0         0         0         0
%                     11.3000    5.4000    2.2000         0   55.8000    5.5000    1.6000    0.8000         0   11.5000    3.1000    2.8000         0         0         0         0         0         0
%                      4.6000   10.1000   11.0000   13.2000    6.5000    4.0000         0         0    1.1000         0         0         0   36.6000   12.8000         0         0         0         0
%                           0         0   13.4000   32.9000    3.7000    7.8000         0         0    8.5000         0         0         0    7.0000   23.0000    3.7000         0         0         0
%                           0         0    3.1000   11.2000         0         0         0         0    1.0000         0         0         0         0   20.4000   39.5000   21.5000         0    6.9000
%                     30.6000         0         0         0   19.8000    8.1000         0         0         0         0         0    7.1000         0         0         0         0   29.1000    5.3000
% 
%                     cpx1       cpx2      cpx3    cpx4       for        fay      spl1       spl2     spl3     opx1       opx2      opx3    fsp1      fsp2      fsp3       qtz       phl       wat      
% cmp_mem = [        5         0         0         0   82.8000    8.2000         0         0         0         0         5         0         0         0         0         0         0         0
%              23.6000    6.4000         0         0   33.6000    4.7000    2.4000    1.2000         0   25.6000    2.6000         0         0         0         0         0         0         0
%               3.3000    2.3000    1.7000         0   75.9000    4.7000         0         0         0    2.1000    6.2000    3.8000         0         0         0         0         0         0
%                    0   17.1000    8.7000   15.8000    2.0000    1.4000         0         0    1.1000         0         0         0   30.5000   23.4000         0         0         0         0
%                    0         0    5.3000   26.3000         0         0         0         0    6.0000         0         0         0         0   29.7000   29.4000    2.5000         0    0.9000
%              38.9000         0         0         0   12.9000   10.2000         0         0         0         0         0    6.5000         0         0         0         0   28.0000    3.5000];

                    %cpx1       cpx2      cpx3         for        fay      spl1       spl2    opx1       opx2     fsp1      fsp2      fsp3       qtz       phl       wat      
% cmp_mem = [        0         0         0   91.84	  8.16         0         0         0         0         0         0         0          0          0          0
              %  27.72	  4.75	    0.00	33.66	  6.83	    3.17	  0.00	     17.43	    6.44   0.00	      0.00	      0.00	    0.00	  0.00	   0.00
              % 1.8000   26.1000    9.4000    6.7000    7.2000    0.8000    0.1000         0         0   45.3000    2.5000         0         0         0         0
              %      0    1.8000   40.0000    0.1000    2.7000         0   12.5000         0         0    7.1000   35.7000    0.1000         0         0         0
              %      0    3.4000   27.1000         0         0         0    1.0000         0         0         0   42.1000   24.5000    1.5000         0    0.5000
              %      4        0         3         0         0         0         3         4         0         0         0        40         0         5    2.9729 ];

              %cpx1       cpx2     cpx3     for       fay      spl1       spl2       opx1       opx2     fsp1      fsp2      fsp3       qtz       phl       wat
cmp_mem = [        0         0         0   91.84	  8.16         0         0         0         0         0         0         0          0          0          0
               27.72	  4.75	    0.00	33.66	  6.83	    3.17	  0.00	     17.43	    6.44   0.00	      0.00	      0.00	    0.00	  0.00	   0.00
              1.8000   26.1000    9.4000    6.7000    7.2000    0.8000    0.1000         0         0   45.3000    2.5000         0         0         0         0
                    0    1.8000   40.0000    0.1000    2.7000         0   12.5000         0         0    7.1000   35.7000    0.1000         0         0         0
                    0    1.1000   28.3000         0         0         0    3.3000         0         0         0   47.4000   19.7000    4.000         0    0.5000
              11.4000         0    9.6000         0         0         0    2.1000    7.5000         0         0         0   66.5000         0         2    0.8000];


%PHS_frc(:,end+1) = 0; 

cmp_mem_init = round(cmp_mem,1);
cmp_mem_init = cmp_mem_init./sum(cmp_mem_init,2)*100;
indmem = logical(cmp_mem_init);
cmp_mem_best = cmp_mem_init;


cmp_oxd_init = cmp_mem_init*cal.mem_oxd/100;
cmp_oxd_best = cmp_oxd_init;


%Set initial component–endmember conditions for fitting melting point parameters.
T0_init = [1850        1605        1120        1044         865         700];  T0_best = T0_init; % T0_init = [1850   1663    1180   807  600] [1850 1600  1120 1000 900 600]
A_init  = [6.5000    5.6000    5.2000    4.2000    1.9000    2.7000];  A_best =  A_init;  %  A_init  = [6.9    5.7   4.5   1.6    3] [7.4  5.2  5.2  3.2  3.2 2.9]
B_init  = [2.6000    3.1000    2.5000    4.1000    1.5000    2.6000];  B_best =  B_init;  %  B_init  = [3.1551    2.8711    2.7025    2.9255    2.8689][3.2  2.9  2.7  3.1  3.1 2.9]
r_init  = [17.8000    7.7000    3.0000    4.7000    7.3000    8.0000]; r_best =  r_init; %   r_init  = [20.5  11.1   7.9   27.5   7.3][30.7  7.5  3.0  5   5 7.4]


% compose initial parameter guess
m0     = [T0_init.';A_init.';B_init.';r_init.';cmp_mem_init(:).*indmem(:);];

% set function to calculate forward model
% m --> dhat
% dhatFunc  = @(model) OxdFromCmpMem(model,MLTp,SOLp,PHS(:,1),cal);
dhatFunc  = @(model) ModelFitNP(model,Tmp,Prs,SYS_oxdp,PHS_frc,cal,[0.5,5,1.75,1e-3]);

% test fit function for initial guess
[~,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmpfit,MLT_cmpfit,SYS_cmpfit] = dhatFunc(m0);

% evaluate melting points as function of pressure
% if isfield(cal,'Tsol'); cal = rmfield(cal,{'Tsol' 'Tliq'}); end
% %PP         = linspace(0.001,max(Psl)*3,50).';
% PP         = linspace(0.001,1.25*3,50).';
% var.m      = ones(size(PP))/2; var.x = var.m; var.f = 0*var.m;
% cal.T0     = T0_init;
% cal.A      = A_init;
% cal.B      = B_init;
% cal.r      = r_init; 
% cal.dTH2O  = dT_init;
% var.c      = repmat(SYS_cmpfit(1,:).*[ones(1,cal.ncmp-1),0]./sum(SYS_cmpfit(1,1:end-1),2),length(PP),1);   % component fractions [wt]
% var.P      = PP/10;
% var.T      = 1000+PP*5e-8;
% var.H2O    = zeros(size(PP));
% cal.H2Osat = var.H2O+0.01;
% [~,cal]    = leappart(var,cal,'T');
% Tm         = cal.Tm;

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
Niter           = 1e5 ;% 1e4;               % 1e6 number of samples to take; how many samples to take during one Monte Carlo run (how many times you do a random shuffle)
anneal.initstep =   5e-4 ;%1e-3;              % Adjust the step size to achieve a reasonable acceptance ratio of 20–30 per cent. The step size determines how far from the current solution guess the next random perturbation is applied
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
% （****）sigma_TSL =  0.1  * ones(size([Tsol(:);Tliq(:)])); % uncertainty of solidus/liquidus Temp
%sigma = [sigma_MLT;sigma_SOL;sigma_MEM;sigma_PHS;sigma_TSL]; % combine all as in data vector
sigma = [sigma_MLT;sigma_SOL;sigma_MEM;sigma_PHS]; % combine all as in data vector

% load calibration data constraints into data vector
%data   = [MLT_oxdp(:);SOL_oxdp(:);SOL_mem(:);PHS_frc(:);Tsol(:);Tliq(:)];
data   = [MLT_oxdp(:);SOL_oxdp(:);SOL_mem(:);PHS_frc(:)];

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
%[dhat,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmpfit,MLT_cmpfit,SYS_cmpfit,Tsolfit,Tliqfit,~] = dhatFunc(bestfit);
[dhat,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmpfit,MLT_cmpfit,SYS_cmpfit] = dhatFunc(bestfit);
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
% var.H2O    = zeros(size(PP));
% cal.H2Osat = var.H2O+0.001;
% [~,cal]    = leappart(var,cal,'T');
% Tm         = cal.Tm;

% uncomment following line to retrieve probability distributions (slow!)
% Nbins = min(500,Niter/20);
% [ppd_mcmc.m, ppd_mcmc.prob] = CalcPDF(mbnds, models(anneal.burnin:end,:), Nbins);


%% *****  visualise best fit calibration  *********************************

%!!!  adjust desired level of detail to plot then run section  !!!
%     level = 1     only simple line plots
%     level = 2     add Harker diagrams and T-X pseudo-sections
%     level = 3     add mineral systems and display best fit parameters

level = 3;
run('MCMC/PlotFit.m')

%% save and display calibration
save('MtAp_noamp10_calibration_2');