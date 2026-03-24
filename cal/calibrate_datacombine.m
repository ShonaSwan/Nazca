
clear all; close all;

%% Load the two processed .mat files (no paper)
%load('MORB_fr_melt.mat');
load('MORB_fr_melt_new.mat');
melt_PHS_frc = PHS_frc;
melt_PHS_oxd = PHS_oxd;
melt_RHO     = RHO;
melt_hasphs  = hasphs;
melt_phs     = phs;
melt_pts     = pts;
melt_Tmp     = Tmp;
melt_Prs     = Prs;
melt_npts    = npts;
melt_noxd    = size(PHS_oxd,3);

%load('MORB_fr_cryst.mat');
load('MORB_fr_cryst_new.mat');
cryst_PHS_frc = PHS_frc;
cryst_PHS_oxd = PHS_oxd;
cryst_RHO     = RHO;
cryst_hasphs  = hasphs;
cryst_phs     = phs;
cryst_pts     = pts;
cryst_Tmp     = Tmp;
cryst_Prs     = Prs;
cryst_npts    = npts;

%% 2. Union of all unique phases (no paper)
all_phs = unique([melt_phs; cryst_phs], 'stable');
nphs = length(all_phs);

% Map old phase indices
melt_iphs = zeros(length(melt_phs),1);
for i = 1:length(melt_phs)
    melt_iphs(i) = find(strcmp(all_phs, melt_phs(i)));
end

cryst_iphs = zeros(length(cryst_phs),1);
for i = 1:length(cryst_phs)
    cryst_iphs(i) = find(strcmp(all_phs, cryst_phs(i)));
end

%% 3. Combine (no CO2 padding needed — use original oxide count)
npts_total = melt_npts + cryst_npts;
noxd = melt_noxd;  

PHS_frc = zeros(npts_total, nphs);
RHO     = zeros(npts_total, nphs);
hasphs  = zeros(npts_total, nphs);
PHS_oxd = zeros(npts_total, nphs, noxd);

% Melting
PHS_frc(1:melt_npts, melt_iphs)                 = melt_PHS_frc;
RHO(1:melt_npts, melt_iphs)                     = melt_RHO;
hasphs(1:melt_npts, melt_iphs)                  = melt_hasphs;
PHS_oxd(1:melt_npts, melt_iphs, :)              = melt_PHS_oxd;

% Crystallisation
PHS_frc(melt_npts+1:melt_npts+cryst_npts, cryst_iphs) = cryst_PHS_frc;
RHO(melt_npts+1:melt_npts+cryst_npts, cryst_iphs)     = cryst_RHO;
hasphs(melt_npts+1:melt_npts+cryst_npts, cryst_iphs)  = cryst_hasphs;
PHS_oxd(melt_npts+1:melt_npts+cryst_npts, cryst_iphs, :) = cryst_PHS_oxd;

%% 4. Update points
pts = [melt_pts; cryst_pts + max(melt_pts)];
Tmp = [melt_Tmp; cryst_Tmp];
Prs = [melt_Prs; cryst_Prs];
npts = npts_total;
phs  = all_phs;

%% 5. Derived calculations
PHS_oxdp = PHS_oxd;

MLT_oxd = squeeze(PHS_oxd(:,1,:));

SOL_oxd = zeros(npts,noxd); wt = zeros(size(SOL_oxd)) + eps;
for iph = 2:nphs
    SOL_oxd = SOL_oxd + squeeze(PHS_oxd(:,iph,:)).*PHS_frc(:,iph);
    wt = wt + PHS_frc(:,iph);
end
SOL_oxd = SOL_oxd ./ wt;

SYS_oxd = zeros(npts,noxd); wt = zeros(size(SYS_oxd)) + eps;
for iph = 1:nphs
    SYS_oxd = SYS_oxd + squeeze(PHS_oxd(:,iph,:)).*PHS_frc(:,iph);
    wt = wt + PHS_frc(:,iph);
end
SYS_oxd = SYS_oxd ./ wt;

%% 6. Save
save('MeltandCryst_new.mat', 'PHS_frc','PHS_oxd','PHS_oxdp','MLT_oxd','SOL_oxd','SYS_oxd',...
     'RHO','phs','hasphs','pts','Tmp','Prs','npts','nphs');

fprintf('COMBINED DATASET SAVED (melting + crystallisation )\n');
fprintf('Total points: %d | Total phases: %d | Oxides: %d\n', npts, nphs, noxd);
fprintf('Ready for PCA / mineral end-member extraction!\n');