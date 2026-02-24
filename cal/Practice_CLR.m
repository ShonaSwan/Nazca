
%% Choose normalisation

%Define some data that looks like oxides values 
X = [
    0.41, 0.12, 0.13, 0.34;
    0.25, 0.30, 0.20, 0.25;
    0.10, 0.50, 0.15, 0.25;
    0.60, 0.05, 0.30, 0.05;
    0.30, 0.40, 0.10, 0.20
];

%% Ask which nornalisation to use ?

scale = upper(input('MSD or CLR? ', 's'));

switch scale
    case 'MSD' % Mean Standard Deviation 
        fprintf('\n=== Row-wise MSD (z-score) ===\n');
        row_mean = mean(X,2);
        row_std  = std(X,0,2);
        Xs = (X - row_mean) ./ row_std;

    case 'CLR' % Centred Log Ratio 
        fprintf('\n=== CLR per row ===\n');
        g = geomean(X,2);
        Xs = log(X ./ g);

    otherwise
        error('Choose MSD or CLR');
end

disp('Normalized Xs:'); disp(Xs);

%% PCA on Xs
fprintf('\n=== PCA on normalized Xs ===\n');
[coeff, score, latent] = pca(Xs, 'Centered', 'on');   % coeff = loadings (DGN.PA), score = scores (DGN.PC)

% Explained variance
explained = latent / sum(latent) * 100;
fprintf('Explained variance by PC: %.1f%%, %.1f%%, ...\n', explained(1), explained(2));

%% Choose number of components to keep (p-1)
n_pc_keep = 2;   % (1,2,3...)


%% Reconstruct in reduced space (simulate Xp)
reduced_score = score(:,1:n_pc_keep);  % truncated scores
reconstructed = reduced_score * coeff(:,1:n_pc_keep)';  % back to normalised space (Xs approx)

% Undo PCA centering (add back column means of Xs)
mu_Xs = mean(Xs,1);
Xs_approx = reconstructed + mu_Xs;

fprintf('Reconstructed in normalized space (Xs_approx):\n');
disp(Xs_approx);

%% Back-transform to original space (FMC)
fprintf('Back-transformed to original space (Xp)');

switch scale
    case 'MSD' % Mean Standard Deviation 
        Xp = Xs_approx .* row_std + row_mean;

    case 'CLR' % Centred Log Ratio 
        Xp = exp(Xs_approx + log(g));
        Xp = Xp ./ sum(Xp, 2);   % force sum-to-1
end

disp('Reconstructed Xp (should be close to original X):');
disp(Xp);

fprintf('Max difference vs original X: %.2e\n', max(abs(X(:) - Xp(:))));
fprintf('Row sums of Xp:'); disp(sum(Xp,2));
