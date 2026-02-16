% ***** TRACE & ISOTOPE GEOCHEMISTRY  *************************************

% *****  Trace Elements  **************************************************

bnd_TRC = zeros(Nz,Nx,cal.ntrc);
Ktrc    = zeros(Nz,Nx,cal.ntrc);
for i = 1:cal.ntrc
    
    % update bulk partitioning coefficients
    for j=1:cal.nmem; Ktrc(:,:,i) = Ktrc(:,:,i) + cal.Ktrc_mem(i,j) .* cx_mem(:,:,j)./100; end

    % update trace element phase compositions
    trcm(:,:,i) = trc(:,:,i)./(m + x.*Ktrc(:,:,i));
    trcx(:,:,i) = trc(:,:,i)./(m./Ktrc(:,:,i) + x);
end

% get trace element advection
[advn_TRCm,qz_advn_TRCm,qx_advn_TRCm] = advect(M.*trcm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);
[advn_TRCx,qz_advn_TRCx,qx_advn_TRCx] = advect(X.*trcx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);

% major component dispersion
[diff_TRCm,qz_diff_TRCm,qx_diff_TRCm] = diffus(trcm,M.*kd,h,[1,2],BCD);
[diff_TRCx,qz_diff_TRCx,qx_diff_TRCx] = diffus(trcx,X.*kx,h,[1,2],BCD);

bnd_TRC = zeros(size(TRC));
if ~isnan(trcwall(1)); bnd_TRC = bnd_TRC + ((trcwall(1,:,:).*rho)-TRC)./(tau_T+dt) .* topshape; end
if ~isnan(trcwall(2)); bnd_TRC = bnd_TRC + ((trcwall(2,:,:).*rho)-TRC)./(tau_T+dt) .* botshape; end

% get total rate of change
dTRCdt = - advn_TRCm - advn_TRCx + diff_TRCm + diff_TRCx + bnd_TRC + Gemt + Gext + Gint;

% residual of trace element evolution
res_TRC = (a1*TRC-a2*TRCo-a3*TRCoo) - (b1*dTRCdt + b2*dTRCdto + b3*dTRCdtoo)*dt;

% semi-implicit update of trace element density
[TRC,GHST.TRC,FHST.TRC,specrad.TRC] = iterate(TRC,res_TRC,specrad.TRC,GHST.TRC,FHST.TRC,itpar,iter);

% convert from densites to concentrations
trc = TRC./RHO;

