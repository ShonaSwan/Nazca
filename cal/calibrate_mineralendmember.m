

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
Fs = {'FontSize',14};
FS = {'FontSize',16};
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


 %Load combined data

 %load Combine_MOR.mat
 load M_Only_MOR.mat 

%liq = 1; ol = 2;  pl = 3;  cpx = 4; spl = 5; bi= 6; q = 7; 
liq = 1; opx = 2; cpx = 3; spl = 4; ol = 5; phl = 6;  pl = 7; 
%liq = 1; opx = 2; cpx = 3; spl = 4; ol = 5; pl = 6; bi= 7; q = 8; phl = 9;  
% !!!  Run Section as is, follow unmix prompts on command line  !!!
cal_MORB_combo;  % read cal.oxdStr from calibration file


DATA.PRJCT  = 'cal';
figno = 100;
%% 

% %***PARAMETER FROM ABOVE SECTION 
hasoxd = logical(squeeze(sum(PHS_oxd,1)));
noxd = 10;
Si = 1; Ti = 2; Al = 3; Cr = 4; Fe = 5; Mg = 6; Ca = 7; Na = 8; K = 9; H = 10;      % set shortcut oxide indices

% initialise lists
PHS_nmem = zeros(nphs,1);
MEM_oxd = [];


% loop through all solid phases
for iph=2:nphs %-1

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
      
% % add quartz as second-last end-member
% MEM_oxd = [MEM_oxd; 100.0, zeros(1,noxd-1)];

% add phlogopite as last end-member
%MEM_oxd = [MEM_oxd; 40.50    2.60   15.90         0      3.00     23.70        0          0    10.700       3.6];

% record final end-member count
% nmem = sum(PHS_nmem) + 1;

% PHS_nmem(1) = nmem;

% !!!  set MEM_oxd => cal.mem_oxd in cal_MORB.m  !!!

%% *****  use end-members to project reduced solid, melt, system compositions
% Compute endmember system oxide composition to approximate the real system. Original values are the true oxide composition

% !!! update calibration file name on following line, then Run Section  !!!
% cal_MORB;  % read cal.mem_oxd from calibration file


% %                SiO2     TiO2      Al2O3     Cr2O3    FeO         MgO       CaO      Na2O     K2O        H2O
%   cal.mem_oxd = [ 30.8700         0         0         0   63.4100    5.7300         0         0         0         0
%                   41.0700         0         0         0    9.3600   49.5700         0         0         0         0
%                   47.9100         0   33.4200         0         0         0   16.4000    2.0500    0.2200         0
%                   64.4300         0   22.2000         0         0         0    3.2800    9.4300    0.6600         0
%                   65.2700         0   20.0900         0         0         0    1.3500    5.8000    7.4900         0
%                   53.0200    0.4900    3.0700         0    4.2500   18.5600   20.2300    0.3800         0         0
%                   50.7700    1.1700    0.6300         0   27.0800    5.8300   12.8800    1.6400         0         0
%                         0         0   41.9800   28.0400   10.3700   19.6100         0         0         0         0
%                         0   27.7400         0         0   71.5500    0.7100         0         0         0         0
%                   35.9500    3.4000   11.6900         0   30.3400    5.9800         0         0    9.7000    2.9400
%                  100.0000         0         0         0         0         0         0         0         0         0];
% 
%  nmem = 11;
%  PHS_nmem = [11 2 3 2 2 1 1];

%                   SiO2     TiO2      Al2O3     Cr2O3        FeO         MgO       CaO      Na2O     K2O        H2O
  cal.mem_oxd = [    54.6300         0    5.6000         0    5.7800   32.1900    1.6600         0         0         0
                     53.5300         0    7.0700         0    6.4300   30.9600    2.0100         0         0         0
                     55.4100         0    3.5400         0    6.1400   31.9400    2.9700         0         0         0
                     
                     53.3600    0.6100    6.5500         0    4.6500   23.9000   10.0600    0.8000         0         0
                     52.9200    0.4900    4.0700         0    4.0500   19.5600   17.2300    0.3800         0         0
                     53.8300    0.1000    2.8000         0    3.7000   19.9900   19.5800         0         0         0
                           0         0   41.9800   28.0400   10.3700   19.6100         0         0         0         0
                     41.0700         0         0         0    9.3600   49.5700         0         0         0         0
                     30.8700         0         0         0   63.4100    5.7300         0         0         0         0
                     40.5000    2.6000   15.9000         0    3.0000   23.7000         0         0   10.7000    3.6000
                     47.7600         0   33.5400         0         0         0   16.5300    2.0400    0.1300         0
                     64.1400         0   22.4100         0         0         0    3.5200    9.3300    0.6000         0];
  nmem = 12;
 PHS_nmem = [12 3 3 1 2 1 2];


% % % %               SiO2    TiO2    Al2O3   Cr2O3    FeO     MgO     CaO    Na2O     K2O     H2O
%   cal.mem_oxd = [  53.5400         0    7.0800         0    6.4200   30.9400    2.0200         0         0         0
%                    54.6600         0    5.6600         0    5.8200   32.1100    1.7500         0         0         0
%                    55.2300         0    3.8500         0    6.1900   31.8900    2.8400         0         0         0
%                    53.3600    0.6100    6.5500         0    4.6500   23.9000   10.0600    0.8000         0         0
%                    53.8300    0.1000    2.8000         0    3.7000   19.9900   19.5800         0         0         0
%                    53.0200    0.4900    3.0700         0    4.2500   18.5600   20.2300    0.3800         0         0
%                    50.7700    1.1700    0.6300         0   27.0800    5.8300   12.8800    1.6400         0         0
%                          0         0   41.9800   28.0400   10.3700   19.6100         0         0         0         0
%                          0   27.7400         0         0   71.5500    0.7100         0         0         0         0
%                    41.0700         0         0         0    9.3600   49.5700         0         0         0         0
%                    30.8700         0         0         0   63.4000    5.7300         0         0         0         0
%                    47.7600         0   33.5400         0         0         0   16.5300    2.0400    0.1300         0
%                    64.1400         0   22.4100         0         0         0    3.5200    9.3300    0.6000         0
%                    66.1100         0   19.7300         0         0         0    0.8500    5.7000    7.6100         0
%                    35.9500    3.4000   11.6900         0   30.3400    5.9800         0         0    9.7000    2.9400
%                   100.0000         0         0         0         0         0         0         0         0         0
%                    39.0000    2.0000   11.0000    1.0000    4.0000   24.0000         0         0   11.0000    8.0000];
%   nmem = 17;
%  PHS_nmem = [17 3 4 2 2 3 1 1 1];

% extract melt phase end-member composition and project back to 
% reduced oxide composition
PHS_mem = zeros(npts,nphs,nmem);
MLT_mem = zeros(npts,nmem);
SOL_mem = zeros(npts,nmem);
MEM_oxd = cal.mem_oxd;
kmem = 1;
for iph = 1:nphs %-1
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


X = [MLT_oxdp;SOL_oxdp]; 
X = X./sum(X,2);

% if more than 2 oxides,
% use unmix tool to perform PCA, end-member extraction   （principal component analysis）
DATA.VNAMES = cal.oxdStr; %%SOMEtime K is empty and seems to cause issues 
DATA.SNAMES = {};
DATA.X      = X;
unmix

ncmp = DGN.p;

MLT_oxdr =  max(0,Xp(   0+(1:npts),:))*100; MLT_oxdr = MLT_oxdr./sum(MLT_oxdr,2)*100;
SOL_oxdr =  max(0,Xp(npts+(1:npts),:))*100; SOL_oxdr = SOL_oxdr./sum(SOL_oxdr,2)*100;
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

% !!! update calibration file name on following line, then Run Section  !!!
% cal_MORB;

% plot selected end-member and projected mineral system compositions
kmem = 1;
for iph=2:nphs-2

    iox = find(hasoxd(iph,:)==1);
    nox = length(iox);

    figure(figno); clf; figno=figno+1;
    set(gcf,'Color','w');

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
set(gcf,'Color','w');
spz = ceil(sqrt(noxd));
spx = ceil((noxd-1)/spz);

kk = 2;

ioxd = [1 2 3 4 5 6 7 8 9 10];
for ix = 1:spx
    for iz = 1:spz
        if kk<= noxd
            subplot(spz,spx,kk-1);
            scatter(MLT_oxd (:,ioxd(1)),MLT_oxd (:,ioxd(kk)),25,Tmp,'o'); colormap('copper'); axis tight; hold on
            scatter(SOL_oxd (:,ioxd(1)),SOL_oxd (:,ioxd(kk)),25,Tmp,'s');colormap('copper'); axis tight; hold on
            scatter(SYS_oxd (:,ioxd(1)),SYS_oxd (:,ioxd(kk)),25,Tmp,'d');
            scatter(MLT_oxdp(:,ioxd(1)),MLT_oxdp(:,ioxd(kk)),25,Tmp,'o','filled');
            scatter(SOL_oxdp(:,ioxd(1)),SOL_oxdp(:,ioxd(kk)),25,Tmp,'s','filled');
            scatter(SYS_oxdp(:,ioxd(1)),SYS_oxdp(:,ioxd(kk)),25,Tmp,'d','filled');
             for icp = 1:ncmp
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
save('MORB_FM_Final');





