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


%% *****  load in data from paper ******************************************
filename = './Schmidt_etal_data.csv';  
uiopen(filename,1)
% *** we do not select liq=100 or liq < 10 wt% ? need to decide limits for fractional melting ? 


%% *****  unpack calibration data  ****************************************

% !!!  update table name on following line, then Run Section  !!!
DAT = Schmidt_etal_data;       % residual fraction = 10%        % table name must correspond to table header above
DAT.Properties.VariableNames{'point___'} = 'point'; DAT.Properties.VariableNames{'T__C_'}    = 'TC';DAT.Properties.VariableNames{'P_kbar_'} = 'Pkbar';DAT.Properties.VariableNames{'mode_wt__'} = 'modewt'; DAT.Properties.VariableNames{'density_kg_m3_'} = 'densitykgm3';

DAT.Properties.VariableNames{'SiO2'}  = 'SiO2wt';DAT.Properties.VariableNames{'TiO2'}  = 'TiO2wt';DAT.Properties.VariableNames{'Al2O3'} = 'Al2O3wt';DAT.Properties.VariableNames{'Cr2O3'} = 'Cr2O3wt';DAT.Properties.VariableNames{'FeO'}   = 'FeOwt'; DAT.Properties.VariableNames{'MgO'}   = 'MgOwt'; DAT.Properties.VariableNames{'CaO'}   = 'CaOwt'; DAT.Properties.VariableNames{'Na2O'}  = 'Na2Owt'; DAT.Properties.VariableNames{'K2O'}   = 'K2Owt';DAT.Properties.VariableNames{'MnO'}   = 'MnOwt';DAT.Properties.VariableNames{'NiO'}   = 'NiOwt';DAT.Properties.VariableNames{'P2O5'}   = 'P2O5wt'; DAT.Properties.VariableNames{'H2O'}   = 'H2Owt'; DAT.Properties.VariableNames{'CO2'}   = 'CO2wt';

% ***** NORMALISE: remove MnO, NiO, P2O5 and renormalize to 100 wt% ********
keep_oxd = {'SiO2wt','TiO2wt','Al2O3wt','Cr2O3wt','FeOwt','MgOwt', ...
           'CaOwt','Na2Owt','K2Owt','H2Owt','CO2wt'};

sum_kept = sum(table2array(DAT(:,keep_oxd)), 2);

for i = 1:height(DAT)
    if sum_kept(i) > 0
        DAT{i,keep_oxd} = DAT{i,keep_oxd} ./ sum_kept(i) * 100;
    end
end

DAT = removevars(DAT, {'MnOwt','NiOwt','P2O5wt'});


%% load phase names in order of appearance, liq first
phs = unique(string(DAT.phase),'stable');                                  % load phase list  
nphs = length(phs);                                                        % record number of phases
iliq = find(strcmp(phs,'liq'));                                            % ensure liq comes first
iphs = 1:nphs; iphs(iliq) = []; iphs = [iliq,iphs];
phs  = phs(iphs);


liq = 1; cpx = 2; g = 3;  ol = 4; opx = 5;           % set shortcut phase indices


% set oxide list in preferred sequence   Reorder the oxides
oxd = ["SiO2";"TiO2";"Al2O3";"Cr2O3";"FeO";"MgO";"CaO";"Na2O";"K2O";"CO2";"H2O"];
noxd = length(oxd);                                                        % record number of oxides
  
unique_points = unique(DAT.point, 'stable');
npts = length(unique_points);   
pts = unique_points;            

Tmp = zeros(npts,1);
Prs = zeros(npts,1);
for ip = 1:npts
    current_point = unique_points(ip);
    point_rows = ismember(DAT.point, current_point);  % safer for numeric
    if any(point_rows)
        first_row = find(point_rows, 1, 'first');
        Tmp(ip) = DAT.TC(first_row);
        Prs(ip) = DAT.Pkbar(first_row);
    else
        warning(['No rows found for point ' num2str(current_point)]);
    end
end

Si = 1; Ti = 2; Al = 3; Cr = 4; Fe = 5; Mg = 6; Ca = 7; Na = 8; K = 9; C = 10; H = 11;  % % set shortcut oxide indices

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
    PHS_oxd(hasphs(:,iph)==1,iph,:) = table2array(DAT(DAT.phase==phs(iph),{'SiO2wt','TiO2wt','Al2O3wt','Cr2O3wt','FeOwt','MgOwt','CaOwt','Na2Owt','K2Owt','CO2wt','H2Owt'}));
    PHS_oxd(hasphs(:,iph)==1,iph,:) = PHS_oxd(hasphs(:,iph)==1,iph,:)./(sum(PHS_oxd(hasphs(:,iph)==1,iph,:),3)+eps)*100;
end

% merge K2O and Na2O
PHS_oxd(:,:,Na) = PHS_oxd(:,:,Na) + PHS_oxd(:,:,K);
PHS_oxd(:,:,K) = [];
noxd = noxd - 1;

% liq = 1; g = 2; cpx = 3; ol = 4; opx = 5; spl = 6;                   % update phase indices
Si = 1; Ti = 2; Al = 3; Cr = 4; Fe = 5; Mg = 6; Ca = 7; Na = 8; C = 9; H = 10; 
oxd  = ["SiO2";"TiO2";"Al2O3";"Cr2O3";"FeO";"MgO";"CaO";"Na2O";"CO2";"H2O"];


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

save('Schmidt_etal_data.mat', 'PHS_frc', 'PHS_oxd', 'PHS_oxdp', 'MLT_oxd','SOL_oxd','SYS_oxd','RHO','phs','hasphs','pts','Tmp','Prs','npts','nphs');

fprintf('Done')
