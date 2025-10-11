% create index array to detect melt presence
hasGem = Gem < -eps;

% find top nonzero index in hasm
melt_iz = ones(1, size(hasGem,2));

% Find the first non-zero row index per column
[row, col] = find(hasGem);
[~, first_idx] = unique(col, 'stable');
melt_iz(col(first_idx)) = row(first_idx);

melt_depth = ZZ(melt_iz);