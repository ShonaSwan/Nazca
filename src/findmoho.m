
moho_depth = NaN(1, Nx);

for ix = 1:Nx
    % Get the density values for each column 
    rhocol = rho(:, ix);
    % Smooth out the values over the columns 
    smthmean = smoothdata(rhocol, 'gaussian', 5);
    % Get the vertical gradient 
    gradrho = gradient(smthmean, h);   
    % Remove any negative values 
    gradrho(gradrho < 0) = 0;   
   
    % Find the max peak
    maxgrad = max(gradrho);
    minHeight = 0.5 * maxgrad;    
    minProm   = 0.2 * maxgrad;

    %find peaks 
    [bound_val, bound_zc] = findpeaks(gradrho, Zc, ...
    'MinPeakHeight', minHeight, ...
    'MinPeakProminence', minProm);

    % Uses shallowestet gradient peak
    if ~isempty(bound_zc)
    moho_depth(ix) = min(bound_zc);
    else
    moho_depth(ix) = NaN;
    end

end 

% %--- plot ---
% 
% figure(15); clf;
% subplot(1,2,1);
% imagesc(1:Nx, Zc, rho); axis ij; hold on
% plot(1:Nx, moho_depth, 'k-', 'LineWidth', 2);
% xlabel('Horizontal index'); ylabel('Depth');
% title('Density field with Moho line');
% 
% subplot(1,2,2);
% plot(moho_depth, 'k-');
% xlabel('Horizontal index'); ylabel('Moho depth');
% title('Moho depth vs horizontal position');
% 
% drawnow;