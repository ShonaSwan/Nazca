

%**************************************************************************
%*****  CALIBRATE MULTI-COMPONENT MELTING MODEL  **************************
%**************************************************************************

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


%% *****  simplify mineral systems and extract end-member compositions  ***
% Principal Component Analysis


% Load combined data
%load MeltandCryst.mat %MORB_fr_cryst.mat
load MeltandCryst_new.mat
liq = 1; g = 2; cpx = 3; ol = 4; opx = 5; spl = 6; pl = 7; 
%liq = 1; g = 2; cpx = 3; ol = 4; opx = 5; spl = 6;  


% !!!  Run Section as is, follow unmix prompts on command line  !!!
cal_FR_combo;  % read cal.oxdStr from calibration file

% prep auxiliary parameters
DATA.PRJCT  = 'cal';
figno = 100;
%% 

% %***PARAMETER FROM ABOVE SECTION 
hasoxd = logical(squeeze(sum(PHS_oxd,1)));
noxd = 10;
Si = 1; Ti = 2; Al = 3; Cr = 4; Fe = 5; Mg = 6; Ca = 7; Na = 8; K = 9;H = 10;      % set shortcut oxide indices

% initialise lists
PHS_nmem = zeros(nphs,1);
MEM_oxd = [];


% loop through all solid phases
for iph=2:nphs

    % extract indices and number of oxides present in phase
    iox = find(hasoxd(iph,:)==1);
    nox = length(iox);


    % load phase compositions into data array for analysis
    X = squeeze(PHS_oxd(hasphs(:,iph)==1,iph,hasoxd(iph,:)));
    X = X./sum(X,2);
    % if more than 2 oxides, 
    % use unmix tool to perform PCA, end-member extraction
    if nox>2 && size(X,1) >= size(X,2)
        DATA.normalisation = 'CLR'; % CLR
        DATA.VNAMES = cal.oxdStr(hasoxd(iph,:));
        DATA.SNAMES = {};
        DATA.X      = X;
        unmix
    % if 2 or less oxides use mean composition as pure-phase end-member
    else
        DGN.p = 1;
    end

    if DGN.p == 1
        FExt = mean(X);
    end


    % process external end-members for phase composition
    EMExt = zeros(DGN.p,noxd);
    EMExt(:,hasoxd(iph,:))  = round(max(0,FExt)./sum(max(0,FExt),2)*100,2);
    EMExt(EMExt==max(EMExt,[],2)) = EMExt(EMExt==max(EMExt,[],2)) + 100 - sum(EMExt,2);

    % sort end-members from primitive to evolved
    [~,is] = sort(EMExt(:,Si)-EMExt(:,Mg)-EMExt(:,Al)+EMExt(:,Na)+EMExt(:,K),'ascend');
    EMExt = EMExt(is,:);

    % add processed end-members to list
    MEM_oxd       = [MEM_oxd;EMExt];
    PHS_nmem(iph) = DGN.p;
end
      
% amph1 = [46  1.7  15.2  0   13   10  11.5  4  0.2  0]
% amph2 = [50  0.25  15   0    9   10  14.5  0.5  0  0]
% MEM_oxd = [MEM_oxd; amph1];
% MEM_oxd = [MEM_oxd; amph2];
% n_amph = 2;   
% PHS_nmem(end+1) = n_amph;   % append a new "phase" slot for amphibole


% add water as last end-member
nmem    = sum(PHS_nmem);
MEM_oxd = [MEM_oxd;zeros(1,noxd-1),100.0];

% record final end-member count and display results
nmem    = sum(PHS_nmem)+1;
PHS_nmem(1) = nmem;
formattedDisplayText(MEM_oxd,'NumericFormat','short')

% !!!  set MEM_oxd => cal.mem_oxd in cal_MORB.m  !!!

%% *****  use end-members to project reduced solid, melt, system compositions
% Compute endmember system oxide composition to approximate the real system. Original values are the true oxide composition

% !!! update calibration file name on following line, then Run Section  !!!
% cal_MORB;  % read cal.mem_oxd from calibration file

% %                SiO2       TiO2      Al2O3     Cr2O3      FeO         MgO       CaO      Na2O     K2O        H2O
cal.mem_oxd = [ 42.3500    0.4100   21.4600    2.8100    5.8000   20.9600    6.2100         0         0         0
                42.2200    0.8900   20.6500    2.0700    7.2600   21.0600    5.8500         0         0         0
                42.1500    1.0000   20.3100    2.1300    7.6500   20.9800    5.7800         0         0         0
                41.1900         0         0         0    7.9600   50.6400    0.2100         0         0         0
                36.4500         0         0         0   32.9500   29.9700    0.6300         0         0         0
                51.9700    0.7000    7.7300    0.2100    4.9700   21.8500   11.5600    0.9800    0.0300         0
                56.8900    0.0600    2.1700    0.0900    4.6100   25.2100    9.4400    1.4100    0.1200         0
                52.7600    0.1400    2.8100    0.4000    3.2600   18.4800   22.1500         0         0         0
                51.9600    1.3200    2.7800    0.2100    7.9300   13.0600   21.8800    0.8200    0.0400         0
                52.3400    0.0700    8.1600    0.6200    6.3400   30.3400    2.0200    0.1100         0         0
                54.7000    0.0400    3.4600    1.2400    5.9600   31.6500    2.9500         0         0         0
                57.5100         0    1.4800    0.0800    5.1800   33.9100    1.4300    0.4100         0         0
                      0    0.0700   51.6800   17.6100    8.7900   21.8400         0         0         0         0
                      0         0   29.5900   40.9600   13.6400   15.8100         0         0         0         0
                      0    1.3600   31.5400   31.3000   33.0800    2.7300         0         0         0         0
                47.2900         0   33.8900         0         0         0   16.9300    1.8900         0         0
                50.5900         0   31.6600         0         0         0   14.3200    3.4200    0.0100         0
                      0         0         0         0         0         0         0         0         0  100.0000]; % water (wat)

 nmem = 18;
 PHS_nmem = [18 3 2 4 3 3 2 1];


% extract melt phase end-member composition and project back to 
% reduced oxide composition
PHS_mem = zeros(npts,nphs,nmem);
MLT_mem = zeros(npts,nmem);
SOL_mem = zeros(npts,nmem);
kmem = 1;
for iph = 1:nphs
    imem  = kmem:kmem+min(nmem,PHS_nmem(iph))-1;
    A     = MEM_oxd(imem,1:end).';
    b     = squeeze(PHS_oxd(:,iph,1:end));
    PHS_mem (:,iph,imem) = lsqregcmp(A,b,[0.01 0 1])*100;

    PHS_oxdp(:,iph,:   ) = squeeze(PHS_mem (:,iph,imem))*MEM_oxd(imem,:)/100 .* max(hasoxd);
    
    if iph==1
        MLT_mem  = squeeze(PHS_mem (:,1,:));
        MLT_oxdp = squeeze(PHS_oxdp(:,1,:));
    else
        SOL_mem(:,imem) = squeeze(PHS_mem (:,iph,imem));
        SOL_mem(:,imem) = SOL_mem(:,imem) .* PHS_frc(:,iph)./(100-PHS_frc(:,1)+eps);
        kmem = kmem+PHS_nmem(iph); 
    end

end
SOL_oxdp = SOL_mem*MEM_oxd/100;           % ‘p’ composition that is projected

% reconstitute projected system oxide composition
SYS_oxdp = zeros(npts,noxd);
wt  = zeros(size(SYS_oxdp)) + eps;
for iph = 1:nphs
    SYS_oxdp = SYS_oxdp + squeeze(PHS_oxdp(:,iph,:)).*PHS_frc(:,iph);
    wt  = wt + PHS_frc(:,iph);
end
SYS_oxdp = SYS_oxdp./wt;  % projected system oxide composition

%% *****  reduce dimensionality by selecting number of pseudo-components **  


X = [MLT_oxdp(:,1:end-1);SOL_oxdp(:,1:end-1)]; 
X = X./sum(X,2);

% if more than 2 oxides,
% use unmix tool to perform PCA, end-member extraction   （principal component analysis）
DATA.VNAMES = cal.oxdStr(1:end-1); %%SOMEtime K is empty and seems to cause issues 
DATA.SNAMES = {};
DATA.X      = X;
unmix

ncmp = DGN.p+1;

MLT_oxdr = [max(0,Xp(   0+(1:npts),:))*100,MLT_oxdp(:,end)]; MLT_oxdr = MLT_oxdr./sum(MLT_oxdr,2)*100;
SOL_oxdr = [max(0,Xp(npts+(1:npts),:))*100,SOL_oxdp(:,end)]; SOL_oxdr = SOL_oxdr./sum(SOL_oxdr,2)*100;
SYS_oxdr = PHS_frc(:,1)/100 .* MLT_oxdr + (1-PHS_frc(:,1)/100) .* SOL_oxdr;

% process internal end-members for phase composition
EMInt = round(max(0,FInt)./sum(max(0,FInt),2)*100,2);
EMInt(EMInt==max(EMInt,[],2)) = EMInt(EMInt==max(EMInt,[],2)) + 100 - sum(EMInt,2);

% sort end-members from primitive to evolved
[~,is] = sort(EMInt(:,Si)-EMInt(:,Mg)+EMInt(:,Na),'ascend');
EMInt  = EMInt(is,:);

% process external end-members for phase composition
EMExt = round(max(0,FExt)./sum(max(0,FExt),2)*100,2);
EMExt(EMExt==max(EMExt,[],2)) = EMExt(EMExt==max(EMExt,[],2)) + 100 - sum(EMExt,2);

% sort end-members from primitive to evolved
[~,is] = sort(EMExt(:,Si)-EMExt(:,Mg)+EMExt(:,Na),'ascend');
EMExt = EMExt(is,:);


%% *****  visualised calibrated end-member, phase compositions  ***********
% These plots show the differences between the original and projected mineral compositions, and are used to judge whether the fit is adequate or if a four endmember model is needed.
%（make a judgment, is this good enough？）
% !!! update calibration file name on following line, then Run Section  !!!
% cal_MORB;

% plot selected end-member and projected mineral system compositions
kmem = 1;
for iph=2:nphs

    iox = find(hasoxd(iph,:)==1);
    nox = length(iox);

    figure(figno); clf; figno=figno+1;

    spz = ceil(sqrt(nox-1));
    spx = ceil((nox-1)/spz);

    kk = 2;
    for ix = 1:spx
        for iz = 1:spz
            if kk<=nox
                subplot(spz,spx,kk-1);
                scatter(squeeze(PHS_oxd (hasphs(:,iph)==1,iph,iox(1))),squeeze(PHS_oxd (hasphs(:,iph)==1,iph,iox(kk))),25,Tmp(hasphs(:,iph)==1)); colormap('copper'); hold on
                scatter(squeeze(PHS_oxdp(hasphs(:,iph)==1,iph,iox(1))),squeeze(PHS_oxdp(hasphs(:,iph)==1,iph,iox(kk))),25,Tmp(hasphs(:,iph)==1),'filled');
                for iem = kmem:kmem+PHS_nmem(iph)-1
                    scatter(MEM_oxd    (iem,iox(1)),MEM_oxd    (iem,iox(kk)),200,'kh','filled');
                    % scatter(cal.mem_oxd(iem,iox(1)),cal.mem_oxd(iem,iox(kk)),200,'kh','filled');
                end
                xlabel(cal.oxdStr(iox(1 )),FS{:},TX{:})
                ylabel(cal.oxdStr(iox(kk)),FS{:},TX{:})
                kk = kk+1;
            else 
                break;
            end
        end
    end
    sgtitle([phs{iph},' projected'],FL{:},TX{:});
    kmem = kmem+PHS_nmem(iph);
    drawnow;
end

% plot fitted liquid, solid, mixture compositions
figure(figno); clf; figno=figno+1;

spz = ceil(sqrt(noxd-1));
spx = ceil((noxd-1)/spz);

kk = 2;

ioxd = [1 2 3 4 5 6 7 8 9];
for ix = 1:spx
    for iz = 1:spz
        if kk<= noxd-1
            subplot(spz,spx,kk-1);
            scatter(MLT_oxd (:,ioxd(1)),MLT_oxd (:,ioxd(kk)),25,Tmp,'o'); colormap('copper'); axis tight; hold on
            scatter(SOL_oxd (:,ioxd(1)),SOL_oxd (:,ioxd(kk)),25,Tmp,'s');colormap('copper'); axis tight; hold on
            scatter(SYS_oxd (:,ioxd(1)),SYS_oxd (:,ioxd(kk)),25,Tmp,'d');
            scatter(MLT_oxdp(:,ioxd(1)),MLT_oxdp(:,ioxd(kk)),25,Tmp,'o','filled');
            scatter(SOL_oxdp(:,ioxd(1)),SOL_oxdp(:,ioxd(kk)),25,Tmp,'s','filled');
            scatter(SYS_oxdp(:,ioxd(1)),SYS_oxdp(:,ioxd(kk)),25,Tmp,'d','filled');
             for icp = 1:ncmp-1
                 if kk<noxd
                     scatter(EMInt(icp,ioxd(1)),EMInt(icp,ioxd(kk)),200,'kh','filled');
                     %scatter(EMExt(icp,ioxd(1)),EMExt(icp,ioxd(kk)),200,'kh');
                 end
             end
            if kk==noxd; legend([{'orig. mlt'},{'orig. sol'},{'orig. sys'},{'proj. sol'},{'proj. mlt'},{'proj. sys'},{'init cmp'}],Fs{:},TX{:},LO{:}); end
            xlabel(cal.oxdStr(ioxd( 1)),FS{:},TX{:})
            ylabel(cal.oxdStr(ioxd(kk)),FS{:},TX{:})
            set(gca,Fs{:},TL{:});
            kk = kk+1;
        else
            break;
        end
    end
end
sgtitle('MLT \& SOL projected',FL{:},TX{:})
drawnow



%% *****  save progress for later use  ************************************

% !!!  Run Section to save calibrated end-members and reduced compositions  !!!
close all;
save('MORB_Comb');


