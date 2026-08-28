% % ***** TRACE & ISOTOPE GEOCHEMISTRY  *************************************
% 
% % *****  Trace Elements  **************************************************
% 
% % major component advection
% [advn_TRCm,qz_advn_TRCm,qx_advn_TRCm] = advect(TRCm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % melt  advection
% [advn_TRCx,qz_advn_TRCx,qx_advn_TRCx] = advect(TRCx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % solid advection
% 
% % major component dispersion
% [diff_TRCm,qz_diff_TRCm,qx_diff_TRCm] = diffus(m.*trcm, rho.*kd, h, [1,2], BCD); 
% [diff_TRCx,qz_diff_TRCx,qx_diff_TRCx] = diffus(x.*trcx, rho.*kx, h, [1,2], BCD);
% 
% % Boundary conditions
% bnd_TRC = zeros(size(TRCx));
% if ~isnan(trcwall(1)); bnd_TRC = bnd_TRC + ((trcwall(1,:,:).*rho)-TRCx)./(tau_T+dt) .* topshape; end
% if ~isnan(trcwall(2)); bnd_TRC = bnd_TRC + ((trcwall(2,:,:).*rho)-TRCx)./(tau_T+dt) .* botshape; end
% 
% % Total rate of change
% dTRCmdt = - advn_TRCm + diff_TRCm + Gmtrc + Gemt + Gint;         
% dTRCxdt = - advn_TRCx + diff_TRCx + Gxtrc + Gext        + bnd_TRC;
% 
% % Time integration
% res_TRCm = (a1*TRCm-a2*TRCmo-a3*TRCmoo) - (b1*dTRCmdt + b2*dTRCmdto + b3*dTRCmdtoo)*dt;
% res_TRCx = (a1*TRCx-a2*TRCxo-a3*TRCxoo) - (b1*dTRCxdt + b2*dTRCxdto + b3*dTRCxdtoo)*dt;
% 
% [trcm,GHST.TRCm,FHST.TRCm,specrad.TRCm] = iterate(trcm,res_TRCm./(M+1e-9.*rho),specrad.TRCm,GHST.TRCm,FHST.TRCm,itpar,iter,[0,nan]);
% [trcx,GHST.TRCx,FHST.TRCx,specrad.TRCx] = iterate(trcx,res_TRCx./(X+1e-9.*rho),specrad.TRCx,GHST.TRCx,FHST.TRCx,itpar,iter,[0,nan]);
% 
% % impose min/max limits on component densities
% trcm = max(0,trcm);
% trcx = max(0,trcx);
% 
% trc  = m.*trcm + x.*trcx;
% TRCm = M.*trcm;
% TRCx = X.*trcx;
% TRC  = TRCm + TRCx;
% 
% % figure(300)
% % subplot(3,1,1)
% % imagesc(Xc,Zc,trcm(:,:,1)); axis equal tight;  box on; cb = colorbar; colormap(ocean)
% % subplot(3,1,2)
% % imagesc(Xc,Zc,trcm(:,:,3)); axis equal tight;  box on; cb = colorbar; colormap(ocean)
% % subplot(3,1,3)
% % imagesc(Xc,Zc,trcm(:,:,4)); axis equal tight;  box on; cb = colorbar; colormap(ocean)
% % drawnow