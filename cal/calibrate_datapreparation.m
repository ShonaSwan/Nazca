

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


%% *****  load calibration data  ******************************************
%load in the MAGEMin data and opens it.(window with table opens, using
%default selection click Import Selection => Import Data) then opens file 
filename = './C_3_T3_ig.csv';  
uiopen(filename,1)
% *** we do not select liq=100 or liq < 10 wt% ? need to decide limits for fractional melting ? 


%% *****  unpack calibration data  ****************************************

% !!!  update table name on following line, then Run Section  !!!
DAT = C_3_T3_ig;       % residual fraction = 10%        % table name must correspond to table header above
DAT.Properties.VariableNames{'point___'} = 'point'; DAT.Properties.VariableNames{'T__C_'}    = 'TC';DAT.Properties.VariableNames{'P_kbar_'} = 'Pkbar';DAT.Properties.VariableNames{'mode_wt__'} = 'modewt'; DAT.Properties.VariableNames{'density_kg_m3_'} = 'densitykgm3';

DAT.Properties.VariableNames{'SiO2_wt__'}  = 'SiO2wt';DAT.Properties.VariableNames{'TiO2_wt__'}  = 'TiO2wt';DAT.Properties.VariableNames{'Al2O3_wt__'} = 'Al2O3wt';DAT.Properties.VariableNames{'Cr2O3_wt__'} = 'Cr2O3wt';DAT.Properties.VariableNames{'FeO_wt__'}   = 'FeOwt'; DAT.Properties.VariableNames{'MgO_wt__'}   = 'MgOwt'; DAT.Properties.VariableNames{'CaO_wt__'}   = 'CaOwt'; DAT.Properties.VariableNames{'Na2O_wt__'}  = 'Na2Owt'; DAT.Properties.VariableNames{'K2O_wt__'}   = 'K2Owt'; DAT.Properties.VariableNames{'H2O_wt__'}   = 'H2Owt';

% load phase names in order of appearance, liq first
phs = unique(string(DAT.phase),'stable');                                  % load phase list
phs(phs=='system') = [];                                                   % discard system
phs(phs=='qfm') = [];                                                      % discard fO2 buffer
phs(phs=='fl') = [];                                                       % discard fluid phase
nphs = length(phs);                                                        % record number of phases
iliq = find(strcmp(phs,'liq'));                                            % ensure liq comes first
iphs = 1:nphs; iphs(iliq) = []; iphs = [iliq,iphs];
phs  = phs(iphs);

liq = 1; cm = 2; olv = 3; opx = 4; spl = 5; cpx = 6; g = 7;              % set shortcut phase indices


% set oxide list in preferred sequence   Reorder the oxides
oxd  = ["SiO2";"TiO2";"Al2O3";"Cr2O3";"FeO";"MgO";"CaO";"Na2O";"K2O";"H2O"];       % set major oxides
noxd = length(oxd);                                                        % record number of oxides

% extract calculation points
pts  = unique(DAT.point,'stable'); offset = min(pts)-1; pts = pts-offset;  % point numbers
Tmp  = unique(DAT.TC,'stable');                                            % point temperatures
Prs  = unique(DAT.Pkbar,'stable');                                         % point pressures
npts = length(pts);                                                        % number of points

Si = 1; Ti = 2; Al = 3; Cr = 4; Fe = 5; Mg = 6; Ca = 7; Na = 8; K = 9; H = 10;      % set shortcut oxide indices

% detect which phases are stable on which points
% Use 0/1 to indicate whether a given mineral phase is present at each point, as not all phases occur at every point. For example, plagioclase and spinel crystallise first, followed by pyroxene.
hasphs = zeros(npts,nphs);
for iph = 1:nphs
    for ipt = 1:npts
        hasphs(ipt,iph) = any(table2array(DAT(DAT.point==ipt+offset,'phase'))==phs(iph));
    end
end

% extract phase fractions in [wt%]
PHS_frc = zeros(npts,nphs);
for iph = 1:nphs
    PHS_frc(hasphs(:,iph)==1,iph) = table2array(DAT(DAT.phase==phs(iph),'modewt'));
end
PHS_frc = PHS_frc./(sum(PHS_frc,2)+eps)*100;

% extract phase densities in [kg/m3]
RHO = zeros(npts,nphs);
for iph = 1:nphs
    RHO(hasphs(:,iph)==1,iph) = table2array(DAT(DAT.phase==phs(iph),'densitykgm3'));
end

% extract phase oxide compositions in [wt%]
PHS_oxd  = zeros(npts,nphs,noxd);
for iph = 1:nphs
    PHS_oxd(hasphs(:,iph)==1,iph,:) = table2array(DAT(DAT.phase==phs(iph),{'SiO2wt','TiO2wt','Al2O3wt','Cr2O3wt','FeOwt','MgOwt','CaOwt','Na2Owt','K2Owt','H2Owt'}));
    PHS_oxd(hasphs(:,iph)==1,iph,:) = PHS_oxd(hasphs(:,iph)==1,iph,:)./(sum(PHS_oxd(hasphs(:,iph)==1,iph,:),3)+eps)*100;
end


% % lump in spl with cm 
PHS_oxd(:,spl,:) = (PHS_frc(:,spl).*PHS_oxd(:,spl,:) + PHS_frc(:,cm ).*PHS_oxd(:,cm ,:)) ./ (PHS_frc(:,spl) + PHS_frc(:,cm ) + eps);
PHS_oxd(:,cm ,:) = [];
RHO(:,spl)       = (PHS_frc(:,spl)+PHS_frc(:,cm))./(PHS_frc(:,spl)./(RHO(:,spl)+eps) + PHS_frc(:,cm)./(RHO(:,cm)+eps) + eps);
RHO(:,cm)        = [];
PHS_frc(:,spl)   =  PHS_frc(:,spl) + PHS_frc(:,cm); 
PHS_frc(:,cm)    = [];
phs(cm)         = [];
hasphs(:,spl)    = max(hasphs(:,spl),hasphs(:,cm));
hasphs(:,cm )    = [];
nphs             = nphs-1;



liq = 1; olv = 2; opx = 3; spl = 4; cpx = 5; g = 6;                   % update phase indices
% *** here spl = spl + cm

% detect which oxides are present in which phases
hasoxd = logical(squeeze(sum(PHS_oxd,1)));

% remove minor oxides from phases (mean<0.50; max<1.0)
for iph = 1:nphs
    ilim = find(mean(squeeze(PHS_oxd(hasphs(:,iph)==1,iph,:)),1)<0.50 & max(squeeze(PHS_oxd(hasphs(:,iph)==1,iph,:)),[],1)<1.00);
    hasoxd(iph,ilim) = false;
    PHS_oxd(:,iph,~hasoxd(iph,:)) = 0;
    PHS_oxd(hasphs(:,iph)==1,iph,:) = PHS_oxd(hasphs(:,iph)==1,iph,:)./(sum(PHS_oxd(hasphs(:,iph)==1,iph,:),3)+eps)*100;
end
PHS_oxdp = PHS_oxd;

% extract melt oxide composition
MLT_oxd = squeeze(PHS_oxd(:,1,:));  % melt oxide composition

% extract solid phase oxide composition
SOL_oxd = zeros(npts,noxd);
wt  = zeros(size(SOL_oxd)) + eps;
for iph = 2:nphs
    SOL_oxd = SOL_oxd + squeeze(PHS_oxd(:,iph,:)).*PHS_frc(:,iph);
    wt  = wt + PHS_frc(:,iph);
end
SOL_oxd = SOL_oxd./wt;  % solid oxide composition

% extract system oxide composition
SYS_oxd = zeros(npts,noxd);
wt  = zeros(size(SYS_oxd)) + eps;
for iph = 1:nphs
    SYS_oxd = SYS_oxd + squeeze(PHS_oxd(:,iph,:)).*PHS_frc(:,iph);
    wt  = wt + PHS_frc(:,iph);
end
SYS_oxd = SYS_oxd./wt;  % system oxide composition

%%  Rename and Save.mat

save('C_3_T3.mat', 'PHS_frc', 'PHS_oxd', 'PHS_oxdp', 'MLT_oxd','SOL_oxd','SYS_oxd','RHO','phs','hasphs','pts','Tmp','Prs','npts','nphs');

fprintf('Done')
