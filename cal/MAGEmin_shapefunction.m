
%Load In CSV file and save it as a table 
Data = readtable('/Users/sas1n24/Documents/Calibration/MORB_eq_melt.csv');

%Find the correct data (relative degree of melting along a vertical)
melts = Data(Data.phase == "liq", :);

%Defining the vertical ? just use the pressure as proxy for depth?
P = melts.P_kbar_;
T = melts.T__C_;
%defining the max depth
Depthmax = max(P);
% normalising so it ranges 0 - 1 (0 = shallow  1= deep)
xi = (P - min(P)) / (max(P) - min(P));

%make a list of all the melt fraction

m = melts.mode_wt__/ 100;  


% define shape function (geometry weighting)
p = 1;                       % 1 = linear triangle, adjust <1 like globe shape (0.5) >1 christmas tree shape 
shapeFunction = xi.^p;         % deepest melt weighted most, shallow melt weighted least

% sum up melt fractions at each point of the MAGEMin model and multiply by
% the width which is proportional to depth 
% i.e., weight incremental melt by geometry

wmelt_contr = m .* shapeFunction; %weighted melt contribution 


% define all the oxide fields from the table 
oxides = {'SiO2_wt__','TiO2_wt__','Al2O3_wt__','Cr2O3_wt__','FeO_wt__','MgO_wt__','CaO_wt__','Na2O_wt__','K2O_wt__','H2O_wt__','O_wt__'};
oxd = zeros(length(P),length(oxides));
for io = 1:length(oxides)
    oxd(:,io)    = melts.(oxides{io});
end

% sum up all these contributions and divide by the sum to get a weighted
% average component 
pooled = zeros(1, numel(oxides));                  
for i = 1:numel(oxides)
    col = oxd(:,i);                     
    pooled(i) = sum(wmelt_contr .* col) / sum(wmelt_contr);
end

exclude = {'Cr2O3_wt__','H2O_wt__','O_wt__'};
keep_idx = ~ismember(oxides, exclude);

oxides_keep = oxides(keep_idx);
pooled_keep = pooled(keep_idx);
pooled_n = pooled_keep ./ sum(pooled_keep) * 100;


%% load in Global Data set from Gale et al., 2014 

opts = detectImportOptions( ...
    '/Users/sas1n24/Documents/Calibration/Gale_et_al_2014.xlsx', ...
    'Sheet', 2);

opts.VariableNamesRange = 'A1';
opts.DataRange          = 'A2';

Gale_et_al = readtable( ...
    '/Users/sas1n24/Documents/Calibration/Gale_et_al_2014.xlsx', ...
    opts);

oxd_Mg_nor = {'Si8_0','Ti8_0','Al8_0','Fe8_0','Ca8_0','Na8_0','K8_0'};
oxd_Fe_nor = {'Si90','Ti90','Al90','Fe90','Mg90','Ca90','Na90','K90'};

oxd_Mg_data = Gale_et_al{1:9563, oxd_Mg_nor};
oxd_Fe_data = Gale_et_al{1:9563, oxd_Fe_nor};
%Plot the data in terms or Harker plots ? then see where our data falls


%% Plotting comparison Global
% Gale Mg data ranges
gale_min = min(oxd_Fe_data,[],1);
gale_max = max(oxd_Fe_data,[],1);
gale_median = median(oxd_Fe_data,1);

%Plot
figure; hold on;

% Patch for range
x_patch = [1:length(oxd_Fe_nor), fliplr(1:length(oxd_Fe_nor))];
y_patch = [gale_min, fliplr(gale_max)];
patch(x_patch, y_patch, 'b', 'FaceAlpha',0.2, 'EdgeColor','none');

% Median line
plot(1:length(oxd_Fe_nor), gale_median, 'b-', 'LineWidth',1.5);

% Your data
scatter(1:length(oxd_Fe_nor), pooled_n, 100, 'r', 'filled');

% Formatting
set(gca, 'XTick', 1:length(oxd_Fe_nor), 'XTickLabel', oxd_Fe_nor);
ylabel('Oxide wt%');
title('Comparison of Fe-normalized Oxides');
grid on;
legend({'Gale data range','Gale median','Your data'}, 'Location','best');

