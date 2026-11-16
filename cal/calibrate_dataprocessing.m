% calibrate_dataprocessing
clear all; close all; clc;

%% 1. Load dataset
load('Frac_M_MOR.mat');

noxd = size(PHS_oxd,3);

%% 2. Ensure phlogopite is included
phs = unique(phs, 'stable');

all_phs = unique(phs, 'stable');

% Add phases only if they do not already exist
new_phases = {'phl'; 'pl'};
all_phs = unique([all_phs(:); new_phases], 'stable');

nphs = length(all_phs);
phs  = all_phs;


%% 3. Expand phase arrays if needed
nphs_old = size(PHS_frc,2);

if nphs > nphs_old
    PHS_frc(:,nphs) = 0;
    RHO(:,nphs) = 0;
    hasphs(:,nphs) = 0;
    PHS_oxd(:,nphs,:) = 0;
end

%% 4. Remove minor oxides
hasoxd = logical(squeeze(sum(PHS_oxd,1)));

for iph = 1:nphs

    present = find(hasphs(:,iph) == 1);

    if isempty(present)
        continue
    end

    % Extract oxide compositions for this phase
    phase_oxd = reshape(PHS_oxd(present,iph,:), ...
                        length(present),noxd);

    % Detect minor oxides
    mean_oxd = mean(phase_oxd,1);
    max_oxd  = max(phase_oxd,[],1);

    ilim = find(mean_oxd < 0.70 & max_oxd < 2.50);

    hasoxd(iph,ilim) = false;

    % Remove minor oxides
    PHS_oxd(:,iph,~hasoxd(iph,:)) = 0;

    % Renormalise oxide compositions to 100 wt%
    phase_oxd = reshape(PHS_oxd(present,iph,:), ...
                        length(present),noxd);

    totals = sum(phase_oxd,2);
    phase_oxd = phase_oxd ./ (totals + eps) * 100;

    PHS_oxd(present,iph,:) = reshape(phase_oxd, ...
                                     length(present),1,noxd);
end

PHS_oxdp = PHS_oxd;

%% 5. Calculate melt oxide composition
% Assumes phase 1 is the melt
MLT_oxd = reshape(PHS_oxd(:,1,:),npts,noxd);

%% 6. Calculate solid oxide composition
SOL_oxd = zeros(npts,noxd);
wt = zeros(npts,1);

for iph = 2:nphs

    phase_oxd = reshape(PHS_oxd(:,iph,:),npts,noxd);

    SOL_oxd = SOL_oxd + phase_oxd .* PHS_frc(:,iph);
    wt = wt + PHS_frc(:,iph);
end

SOL_oxd = SOL_oxd ./ (wt + eps);

%% 7. Calculate system oxide composition
SYS_oxd = zeros(npts,noxd);
wt = zeros(npts,1);

for iph = 1:nphs

    phase_oxd = reshape(PHS_oxd(:,iph,:),npts,noxd);

    SYS_oxd = SYS_oxd + phase_oxd .* PHS_frc(:,iph);
    wt = wt + PHS_frc(:,iph);
end

SYS_oxd = SYS_oxd ./ (wt + eps);

%% 8. Save processed dataset
save('M_Only_MOR.mat', ...
    'PHS_frc','PHS_oxd','PHS_oxdp', ...
    'MLT_oxd','SOL_oxd','SYS_oxd', ...
    'RHO','phs','hasphs','pts','Tmp','Prs','npts','nphs');

%% 9. Display summary
fprintf('M-ONLY DATASET SAVED\n');
fprintf('Total points: %d | Total phases: %d | Oxides: %d\n', ...
        npts,nphs,noxd);
fprintf('Ready for PCA / mineral end-member extraction!\n');