%% load in Global Data set from Gale et al., 2014 

opts = detectImportOptions( ...
    '/Users/sas1n24/Documents/Calibration/Gale_et_al_2014.xlsx', ...
    'Sheet', 2);

opts.VariableNamesRange = 'A1';
opts.DataRange          = 'A2';

Gale_et_al = readtable( ...
    '/Users/sas1n24/Documents/Calibration/Gale_et_al_2014.xlsx', ...
    opts);

% Oxide names within both data sets 
oxd_Mg_nor = {'Si8_0','Ti8_0','Al8_0','Fe8_0','Ca8_0','Na8_0','K8_0'};
oxd_Fe_nor = {'Si90','Ti90','Al90','Fe90','Mg90','Ca90','Na90','K90'};

% Global Dataset Mg 8 normalised
oxd_Mg_data = Gale_et_al{1:9563, oxd_Mg_nor};
% Global Dataset Fe 90 normalised
oxd_Fe_data = Gale_et_al{1:9563, oxd_Fe_nor};
% Galapagos Dataset 
Gala_data = Gale_et_al{5265:5678, oxd_Fe_nor};

%% Load In my calibration data 
opts = detectImportOptions( ...
    '/Users/sas1n24/Documents/Calibration/Calibration_Data_Shapefunction.xlsx', ...
    'Sheet', 1);

opts.VariableNamesRange = 'A3';
opts.DataRange          = 'A4';


Cal_data = readtable( ...
    '/Users/sas1n24/Documents/Calibration/Calibration_Data_Shapefunction.xlsx', ...
    opts);

oxd = {'SiO2','TiO2','Al2O3','FeO','MgO','CaO','Na2O','K2O','H2O','O','Cr2O3'};
oxd_all = Cal_data{1:12, oxd};

%Normalise to match oxides within Fe 90 dataset 
oxd_val = oxd_all(:, 1:8) ./ sum(oxd_all(:,1:8), 2) * 100;

oxd_val_mg = oxd_val;
oxd_val_mg(:,5) = [];



%% Comparison with Global data qualitative 

oxide_names = {'TiO2', 'Al2O3', 'FeO', 'MgO', 'CaO', 'Na2O', 'K2O'}; 

% Loop through oxides
for j = 1:(size(oxd_val, 2) - 1)
    i = j + 1; 
    
    figure;
    hold on;
    x = oxd_Fe_data(:,1);  % SiO2
    y = oxd_Fe_data(:,i);  % Current oxide
    valid_idx = ~isnan(x) & ~isnan(y);  % Remove NaNs if any
    x = x(valid_idx);
    y = y(valid_idx);
    
    k = convhull(x, y);  % Compute convex hull indices
    fill(x(k), y(k), [0.7 0.7 1], ...  % Lighter blue fill
    'FaceAlpha', 0.3, ...  % Semi-transparent
    'EdgeColor', 'b', ...
    'LineWidth', 1.5);
    
    % Plot the average of Gale et al. 
    x_mean = mean(x); 
    y_mean = mean(y);
    scatter(x_mean, y_mean, 150, 'd','MarkerFaceColor', 'b');

    %Use different colours for my points 
    num_points = size(oxd_val, 1);
    colors = lines(num_points); 
    for p = 1:num_points
        scatter(oxd_val(:,1), oxd_val(:,i), 50, ceil((1:size(oxd_val,1))'/8), 'filled')
        colormap(lines(3))
    end
    
    xlabel('SiO2', 'Interpreter', 'none');
    ylabel(oxide_names{j}, 'Interpreter', 'none');
    axis xy tight
    grid on
    legend([{'Gale et al. Spread'}, 'My calibrations'], 'Location', 'best'); 
    
    title(sprintf('Global Harker Plot for %s', oxide_names{j}));
    hold off;
end

%% Comparison with Galapagos data qualitative 

for j = 1:(size(oxd_val, 2) - 1)
    i = j + 1; 
    
    figure;
    hold on;
    x = Gala_data(:,1);  % SiO2
    y = Gala_data(:,i);  % Current oxide
    valid_idx = ~isnan(x) & ~isnan(y);  % Remove NaNs if any
    x = x(valid_idx);
    y = y(valid_idx);
    
    k = convhull(x, y);  % Compute convex hull indices
    fill(x(k), y(k), [1 1 0], ...        % yellow fill
    'FaceAlpha', 0.25, ...
    'EdgeColor', [0.85 0.65 0], ... % darker yellow/orange edge
    'LineWidth', 1.5);
    
    % Plot the average of Gale et al. 
    x_mean = mean(x); 
    y_mean = mean(y);
   scatter(x_mean, y_mean, 150, 'd', ...
    'MarkerFaceColor', [0.9 0.7 0], ...
    'MarkerEdgeColor', 'k');

    %Use different colours for my points 
    num_points = size(oxd_val, 1);
    colors = lines(num_points); 
    for p = 1:num_points
        scatter(oxd_val(:,1), oxd_val(:,i), 50, ceil((1:size(oxd_val,1))'/6), 'filled')
        colormap(lines(3))
    end
    
    xlabel('SiO2', 'Interpreter', 'none');
    ylabel(oxide_names{j}, 'Interpreter', 'none');
    axis xy tight
    grid on
    legend([{'Gale et al. Galapagos'}, 'My calibrations'], 'Location', 'best'); 
    
    title(sprintf('Galapagos Harker Plot for %s', oxide_names{j}));
    hold off;
end


%% Comparison with Global data quantitative

eps_val = 1e-6; % small enough to not affect compositions
oxd_Fe_data(oxd_Fe_data <= 0) = eps_val;
oxd_val(oxd_val <= 0) = eps_val;
idx = 1:size(oxd_Fe_data,2);

% CLR transform global dataset
global_clr = clr_transform(oxd_Fe_data(:,idx));

% Global mean in CLR space
mu_clr = mean(global_clr, 1, 'omitnan');

num_cal = size(oxd_val,1);
score = zeros(num_cal,1);

% CLR transform candidate compositions
val_clr = clr_transform(oxd_val(:,idx));

for i = 1:num_cal
    d = val_clr(i,:) - mu_clr;
    score(i) = sqrt(sum(d.^2, 'omitnan'));   % Euclidean distance in CLR space
end

[score_sorted, order] = sort(score);
results = table(order, score_sorted, ...
    'VariableNames', {'Calibration_ID','CLR_Distance'});

disp(results)

%% Comparison with Galapagos data quantitative

eps_val = 1e-6; % small enough to not affect compositions
Gala_data(Gala_data <= 0) = eps_val;
oxd_val(oxd_val <= 0) = eps_val;
idx = 1:size(Gala_data,2);

% CLR transform global dataset
global_clr = clr_transform(Gala_data(:,idx));

% Global mean in CLR space
mu_clr = mean(global_clr, 1, 'omitnan');

num_cal = size(oxd_val,1);
score = zeros(num_cal,1);

% CLR transform candidate compositions
val_clr = clr_transform(oxd_val(:,idx));

for i = 1:num_cal
    d = val_clr(i,:) - mu_clr;
    score(i) = sqrt(sum(d.^2, 'omitnan'));   % Euclidean distance in CLR space
end

[score_sorted, order] = sort(score);
results = table(order, score_sorted, ...
    'VariableNames', {'Calibration_ID','CLR_Distance'});

disp(results)


function clr = clr_transform(X)
    % X: samples x oxides (wt%)
    % Normalize to constant sum
    X = X ./ sum(X,2) * 100;

    % Geometric mean per sample
    g = exp(mean(log(X), 2, 'omitnan'));

    % CLR transform
    clr = log(X ./ g);
end

%% Comparison with Global data quantitative FC

eps_val = 1e-6; % small enough to not affect compositions
oxd_Mg_data(oxd_Mg_data <= 0) = eps_val;
oxd_val_mg(oxd_val_mg <= 0) = eps_val;
idx = 1:size(oxd_Mg_data,2);

% CLR transform global dataset
global_clr = clr_transform(oxd_Mg_data(:,idx));

% Global mean in CLR space
mu_clr = mean(global_clr, 1, 'omitnan');

num_cal = size(oxd_val_mg,1);
score = zeros(num_cal,1);

% CLR transform candidate compositions
val_clr = clr_transform(oxd_val_mg(:,idx));

for i = 1:num_cal
    d = val_clr(i,:) - mu_clr;
    score(i) = sqrt(sum(d.^2, 'omitnan'));   % Euclidean distance in CLR space
end

[score_sorted, order] = sort(score);
results = table(order, score_sorted, ...
    'VariableNames', {'Calibration_ID','CLR_Distance'});

disp(results)

