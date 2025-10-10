
moho_depth = NaN(1, Nx);
%Defining composition ratio for moho
%high in crust low in mantle 
c_moho = c_oxd(:,:,1) ./ c_oxd(:,:,5);


for ix = 1:Nx

    % Get the composition values for each column 
    compcol = c_moho(:, ix);
    % Smooth out the values over the columns 
    %smthmean = smoothdata(compcol, 'gaussian', 5);
    % Get the vertical gradient 

    gradcomp = abs(gradient(compcol, h));
    %gradcomp = gradient(smthmean, h);   

    % Remove any negative values 
    gradcomp(gradcomp < 0) = 0;   
   
    % Find the max peak
   
    maxgrad = max(gradcomp);
    minHeight = 0.1 * maxgrad;    %0.5
    minProm   = 0.05 * maxgrad;    %0.2

    %find peaks 
    %[bound_val, bound_zc] = islocalmax(gradcomp, Zc, ...
    [bound_val, bound_zc] = findpeaks(gradcomp, Zc, ...
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

figure(15); clf;
subplot(1,2,1);
imagesc(1:Nx, Zc, c_moho); axis ij; hold on
plot(1:Nx, moho_depth, 'k-', 'LineWidth', 2);
xlabel('Horizontal index'); ylabel('Depth');
title('SiO₂/MgO field with Moho line');
colorbar;
subplot(1,2,2);
plot(moho_depth, 'k-');
xlabel('Horizontal index'); ylabel('Moho depth');
title('Moho depth vs horizontal position');

drawnow;