
horzmean = mean(rho, 2);                          
smthmean = smoothdata(horzmean, 'gaussian', 5);  

gradrho = gradient(smthmean, h);                  
gradrho(gradrho < 0) = 0;                        

maxgrad = max(gradrho);
minHeight = 0.5 * maxgrad;    
minProm   = 0.2 * maxgrad;    

%find peaks 
[bound_val, bound_zc] = findpeaks(gradrho, Zc, ...
    'MinPeakHeight', minHeight, ...
    'MinPeakProminence', minProm);

%pick which peak to use (where the moho is)
if ~isempty(bound_zc)
    moho_depth = min(bound_zc);
else
    moho_depth = NaN;
end

% % --- plot ---
% figure(14); clf;
% sgtitle(['time = ', num2str(time/hr,3), ' [hr]'], TX{:}, FS{:}, 'Color', 'k');
% 
% subplot(1,2,1);
% plot(horzmean, zc, 'b-'); axis ij tight; hold on
% plot(smthmean, zc, 'r-');
% xlabel('\rho'); ylabel('Depth');
% legend('Mean \rho', 'Smoothed');
% 
% subplot(1,2,2);
% plot(gradrho, zc, 'k-'); axis ij tight; hold on
% plot(bound_val, bound_zc, 'ro', 'MarkerFaceColor', 'r');
% xlabel('d\rho/dz'); ylabel('Depth');
% legend('Gradient', 'Detected peaks');
% 
% drawnow;