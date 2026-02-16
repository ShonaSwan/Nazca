%Load In CSV file and save it as a table 
Data = readtable('/Users/sas1n24/Documents/Calibration/C_fc_p_7_ig.csv');

%Find the correct data (relative degree of melting along a vertical)
melts = Data(Data.phase == "liq", :);

%Defining the vertical ? just use the pressure as proxy for depth?
P = melts.P_kbar_;
T = melts.T__C_;
%defining the max depth
Depthmax = max(P);

%make a list of all the melt fraction

m = melts.mode_wt__/ 100;  


% define all the oxide fields from the table 
oxides = {'SiO2_wt__','TiO2_wt__','Al2O3_wt__','Cr2O3_wt__','FeO_wt__','MgO_wt__','CaO_wt__','Na2O_wt__','K2O_wt__','H2O_wt__','O_wt__'};
oxd = zeros(length(P),length(oxides));
for io = 1:length(oxides)
    oxd(:,io)    = melts.(oxides{io});
end

exclude = {'Cr2O3_wt__','H2O_wt__','O_wt__'};
keep_idx = ~ismember(oxides, exclude);

oxides_keep = oxides(keep_idx);
oxd_keep = oxd(:,keep_idx);
oxd_norm = oxd_keep ./ sum(oxd_keep, 2) * 100;

plot_oxd = oxd_norm;
plot_oxd(:,5) = [];

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



%% Comparison with Global data qualitative 

% Oxide names (skip SiO2)
oxide_names = {'TiO2', 'Al2O3', 'FeO', 'CaO', 'Na2O', 'K2O'}; 

% Loop through oxides
for j = 1:size(plot_oxd, 2)-1
    
    figure; hold on;

    % --- Gale background ---
    xg = oxd_Mg_data(:,1);           % SiO2
    yg = oxd_Mg_data(:,j+1);         % current oxide
    valid_idx = ~isnan(xg) & ~isnan(yg);
    xg = xg(valid_idx);
    yg = yg(valid_idx);

    % Convex hull
    k = convhull(xg, yg);
    fill(xg(k), yg(k), [0.7 0.7 1], 'FaceAlpha', 0.3, 'EdgeColor','none');

    % Gale centroid
    scatter(mean(xg), mean(yg), 120, 'b', 'd', 'filled');

    % --- My data ---
    x = plot_oxd(:,1);          % SiO2
    y = plot_oxd(:,j+1);        % current oxide
    in = inpolygon(x, y, xg(k), yg(k));  % logical mask for inside points

    % Keep only points inside the envelope
    x_in = x(in);
    y_in = y(in);

    % Plot only inside points
    scatter(x_in, y_in, 40, 'r', 'filled');

    % Connect inside points with a line
    plot(x_in, y_in, 'r-', 'LineWidth', 2);

    % --- Formatting ---
    xlabel('SiO_2 (wt%)');
    ylabel(oxide_names{j});
    title(['Harker plot: ' oxide_names{j}]);
    grid on;
    axis tight;

    legend({'Gale envelope','Gale centroid','Model inside envelope','Model trend'}, 'Location','best');

    hold off;
end