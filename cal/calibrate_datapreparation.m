

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
filename = './MORB_fr_melt.csv';  
uiopen(filename,1)
% *** we do not select liq=100 or liq < 10 wt% ? need to decide limits for fractional melting ? 


%% *****  unpack calibration data  ****************************************

% !!!  update table name on following line, then Run Section  !!!
DAT = MORB_fr_melt;       % residual fraction = 10%        % table name must correspond to table header above
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
  
%liq = 1; ol = 2; spl = 3; pl = 4; cpx = 5;
liq = 1; g = 2 ;cpx = 3; ol = 4; opx = 5; spl = 6; cm = 7; 

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


% lump in spl with cm 
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

% % merge K2O and Na2O
% PHS_oxd(:,:,Na) = PHS_oxd(:,:,Na) + PHS_oxd(:,:,K);
% PHS_oxd(:,:,K) = [];
% noxd = noxd - 1;

liq = 1; g = 2 ;cpx = 3; ol = 4; opx = 5; spl = 6; 

% detect which oxides are present in which phases
hasoxd = logical(squeeze(sum(PHS_oxd,1)));

% remove minor oxides from phases (mean<0.50; max<1.0)
% for iph = 1:nphs
%     ilim = find(mean(squeeze(PHS_oxd(hasphs(:,iph)==1,iph,:)),1)<0.50 & max(squeeze(PHS_oxd(hasphs(:,iph)==1,iph,:)),[],1)<1.00);
%     hasoxd(iph,ilim) = false;
%     PHS_oxd(:,iph,~hasoxd(iph,:)) = 0;
%     PHS_oxd(hasphs(:,iph)==1,iph,:) = PHS_oxd(hasphs(:,iph)==1,iph,:)./(sum(PHS_oxd(hasphs(:,iph)==1,iph,:),3)+eps)*100;
% end

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

 %% Take the Melt phase and run it through the shape function 
% Melt_fractions = PHS_frc(:,1);                    % melt mode wt% at each point
% Melt_oxd       = squeeze(PHS_oxd(:,liq,:));       % all oxide compositions of the melt at each point
% 
% %%%%% subtract connectivity threshold (1%) and zero out negatives %%%%%%%%
% threshold       = 1;
% extracted_mass  = max(Melt_fractions - threshold, 0);  
% 
% %%%% build the triangular shape function %%%%
% P = zeros(npts,1);
% for i = 1:npts
%     P(i) = DAT.Pkbar(find(DAT.point == i + offset, 1));
% end
% xi            = (P - min(P)) / (max(P) - min(P));   % 0 = shallow, 1 = deepest
% p             = 1;                                  % 1 = simple triangle
% shapeFunction = xi.^p;
% 
% % 3. multiply the extracted mass by the shape and by the oxide values
% wmelt_contr = (extracted_mass / 100) .* shapeFunction;   % convert to fraction
% 
% %%%%%%% sum all points together %%%%%%%
% pooled = zeros(1, noxd);
% for io = 1:noxd
%     pooled(io) = sum(wmelt_contr .* Melt_oxd(:,io)) / (sum(wmelt_contr) + eps);
% end
% pooled_n = pooled ./ sum(pooled) * 100;               % final parental melt composition


%%  Rename and Save.mat

save('MORB_fr_melt_minor.mat', 'PHS_frc', 'PHS_oxd', 'PHS_oxdp', 'MLT_oxd','SOL_oxd','SYS_oxd','RHO','phs','hasphs','pts','Tmp','Prs','npts','nphs');

fprintf('Done')
