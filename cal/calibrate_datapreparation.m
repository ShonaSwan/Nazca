

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

filename = './MOR_Cryst_ig.csv';  
uiopen(filename,1)
% *** we do not select liq=100 or liq < 10 wt% ? need to decide limits for fractional melting ? 


%% *****  unpack calibration data  ****************************************

% !!!  update table name on following line, then Run Section  !!!
DAT = MOR_Cryst_ig;       % residual fraction = 10%        % table name must correspond to table header above
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
  
liq = 1; ol = 2; spl = 3; cpx = 4; pl = 5; ilm = 6; usp = 7; mgt = 8; 
%liq = 1; opx = 2; cpx = 3; ol = 4;  g = 5 ; spl = 6; 

% set oxide list in preferred sequence   Reorder the oxides
oxd  = ["SiO2";"TiO2";"Al2O3";"Cr2O3";"FeO";"MgO";"CaO";"Na2O";"K2O";"H2O"];       % set major oxides
noxd = length(oxd);                                                        % record number of oxides

% extract calculation points
pts  = unique(DAT.point,'stable');   % point numbers
Tmp  = unique(DAT.TC,'stable');                                            % point temperatures
Prs  = unique(DAT.Pkbar,'stable');                                         % point pressures
npts = length(pts);                                                        % number of points

Si = 1; Ti = 2; Al = 3; Cr = 4; Fe = 5; Mg = 6; Ca = 7; Na = 8; K = 9; H = 10;      % set shortcut oxide indices

% detect which phases are stable on which points
% Use 0/1 to indicate whether a given mineral phase is present at each point, as not all phases occur at every point. For example, plagioclase and spinel crystallise first, followed by pyroxene.
hasphs = zeros(npts,nphs);
for iph = 1:nphs
    for ipt = 1:npts
        hasphs(ipt,iph) = any(table2array(DAT(DAT.point==ipt,'phase'))==phs(iph));
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

% lump in mgt with usp and ilm 
PHS_oxd(:,ilm,:) = (PHS_frc(:,ilm).*PHS_oxd(:,ilm,:) + PHS_frc(:,usp).*PHS_oxd(:,usp,:) + PHS_frc(:,mgt).*PHS_oxd(:,mgt,:))./ (PHS_frc(:,ilm) + PHS_frc(:,usp) + PHS_frc(:,mgt) + eps);
PHS_oxd(:, [usp, mgt], :) = [];
RHO(:,ilm) = (PHS_frc(:,ilm) + PHS_frc(:,usp) + PHS_frc(:,mgt))./ (PHS_frc(:,ilm)./(RHO(:,ilm) + eps)+ PHS_frc(:,usp)./(RHO(:,usp) + eps)+ PHS_frc(:,mgt)./(RHO(:,mgt) + eps)+ eps);
RHO(:, [usp, mgt])        = [];
PHS_frc(:,ilm) = PHS_frc(:,ilm) + PHS_frc(:,usp) + PHS_frc(:,mgt);
PHS_frc(:, [usp, mgt])    = [];
phs([usp, mgt])           = [];
hasphs(:,ilm) = max(hasphs(:,ilm), max(hasphs(:,usp), hasphs(:,mgt)));
hasphs(:, [usp, mgt])     = [];
nphs             = nphs-2;


%% Take the Melt phase and run it through the shape function 
Melt_fractions = PHS_frc(:,1);                    % melt mode wt% at each point
Melt_oxd       = squeeze(PHS_oxd(:,liq,:));       % all oxide compositions of the melt at each point

%%%%% subtract connectivity threshold (1%) and zero out negatives %%%%%%%%
threshold       = 1;
extracted_mass  = max(Melt_fractions - threshold, 0);  

%%%% build the triangular shape function %%%%
P = zeros(npts,1);
for i = 1:npts
    P(i) = DAT.Pkbar(find(DAT.point == i, 1));
end
pmax = max(P(extracted_mass>0));
xi            = (P - (2.8)) / (pmax - (2.8));   % 0 = shallow, 1 = deepest
xi(P>pmax)= 0;
p             = 0.3;                                  % 1 = simple triangle
shapeFunction = xi.^p;

% multiply the extracted mass by the shape and by the oxide values
wmelt_contr = (extracted_mass / 100) .* shapeFunction;   % convert to fraction

%%%%%%% sum all points together %%%%%%%
pooled = zeros(1, noxd);
for io = 1:noxd
    pooled(io) = sum(wmelt_contr .* Melt_oxd(:,io)) / (sum(wmelt_contr) + eps);
end
pooled_n = pooled ./ sum(pooled) * 100;               % final parental melt composition
Tot_mel = sum(wmelt_contr .* (1 + cumsum(extracted_mass))) / (sum(wmelt_contr) + eps);


%%  Rename and Save.mat

save('MORB_fr_Cryst.mat', 'PHS_frc', 'PHS_oxd','RHO','phs','hasphs','pts','Tmp','Prs','npts','nphs');

fprintf('Done');
