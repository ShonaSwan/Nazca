
% get vertical gradient of bulk SiO2
[~,dSiO2dz] = gradient(c_oxd(:,:,1));
dSiO2dz(dSiO2dz > 0) = 0;
iz = floor(Nz/4);
dSiO2dz(iz+1:end,:) = repmat(dSiO2dz(iz,:),Nz-iz,1);

% maxgrad = max(-dSiO2dz(:));
% minProm = 0.1 * maxgrad;    
% 
% % find local max along vertical columns 
% bound_zc = islocalmax(-dSiO2dz,1,'MinProminence',minProm,'MaxNumExtrema',1);
% bound_zc(1,~any(bound_zc,1)) = 1;

[ival,imax] = max(-dSiO2dz,[],1);

% retrieve moho location
moho_depth = Zc(imax);

figure(200);
imagesc(Xc,Zc,-dSiO2dz); axis equal tight; colorbar; colormap(ocean); hold on
plot(Xc,moho_depth,'w','LineWidth',1.5)
