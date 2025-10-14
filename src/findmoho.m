%get vertical gradient of bulk SiO2
SiO2ByMgO = c_oxd(:,:,cal.Si)./c_oxd(:,:,cal.Mg);

[~,gradSiO2ByMgO] = gradient(SiO2ByMgO);
gradSiO2ByMgO(gradSiO2ByMgO > 0) = 0;

indLAB    = (T-273.15) < 1200;
[LAB_depth,LAB_iz] = max(ZZ.*indLAB,[],1);

for ix = 1:Nx
gradSiO2ByMgO(LAB_iz(ix):end,ix) = gradSiO2ByMgO(LAB_iz(ix),ix);
end

[~,MOHO_iz] = max(-gradSiO2ByMgO,[],1);

% retrieve moho location
MOHO_depth = Zc(MOHO_iz);

% figure(200);
% imagesc(Xc,Zc,-gradSiO2ByMgO); axis equal tight; colorbar;colormap(colmap); hold on
% plot(Xc,MOHO_depth,'w','LineWidth',1.5)
% plot(Xc,LAB_depth,'w--','LineWidth',1.5)


% moho_depth = NaN(1, Nx);
% 
% %Define compositional ratio of Moho detection
% %high in crust low in mantle 
% c_moho = c_oxd(:,:,1) ./ c_oxd(:,:,5);
% 
% 
% % Get the vertical gradient 
% 
% gradcomp = abs(gradient(c_moho, h, 1));
% 
% 
% % Remove any negative values 
% gradcomp(~isfinite(gradcomp)) = 0;   
% 
% ismax = islocalmax(gradcomp, 1);
% 
% maxgrad_col = max(gradcomp, [], 1, 'omitnan');
% minHeight_col = 0.1 * maxgrad_col;
% minProm_col   = 0.05 * maxgrad_col;
% 
% % --- Apply thresholds ---
% for ix = 1:Nx
%     col = gradcomp(:, ix);
%     peaks = ismax(:, ix) & (col >= minHeight_col(ix));
% 
%     % Find depths corresponding to valid peaks
%     if any(peaks)
%         moho_depth(ix) = Zc(find(peaks, 1, 'first'));  % shallowest significant peak
%     else
%         moho_depth(ix) = NaN;
%     end
% end
% 
% % --- Smooth the detected Moho depth ---
% moho_depth_smooth = smoothdata(moho_depth, 'gaussian', 10, 'omitnan');
% 
% % --- Plot ---
% figure(15); clf;
% 
% % 1️⃣ Field plot with smoothed Moho line
% subplot(1,2,1);
% imagesc(1:Nx, Zc, c_moho); axis ij; hold on
% plot(1:Nx, moho_depth_smooth, 'k-', 'LineWidth', 2);
% xlabel('Horizontal index'); ylabel('Depth');
% title('SiO₂/MgO field with smoothed Moho line');
% colorbar;
% 
% % (optional) show raw line for comparison (in gray dashed)
% plot(1:Nx, moho_depth, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1, 'LineStyle', '--');
% 
% % 2️⃣ Moho depth vs horizontal position
% subplot(1,2,2);
% plot(moho_depth, 'Color', [0.6 0.6 0.6], 'LineStyle', '--'); hold on
% plot(moho_depth_smooth, 'k-', 'LineWidth', 2);
% xlabel('Horizontal index'); ylabel('Moho depth');
% legend('Raw', 'Smoothed');
% title('Moho depth vs horizontal position');
% 
% drawnow;