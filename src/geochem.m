% ***** TRACE & ISOTOPE GEOCHEMISTRY  *************************************

% *****  Trace Elements  **************************************************

bnd_TRC = zeros(Nz,Nx,cal.ntrc);

% Advection + diffusion + disequilibrium compositions
[advn_TRCm,~,~] = advect(M.*trcm, Um(2:end-1,:), Wm(:,2:end-1), h, {ADVN,''}, [1,2], BCA);
[advn_TRCx,~,~] = advect(X.*trcx, Ux(2:end-1,:), Wx(:,2:end-1), h, {ADVN,''}, [1,2], BCA);

[diff_TRCm,~,~] = diffus(m.*trcm, rho.*kd, h, [1,2], BCD);
[diff_TRCx,~,~] = diffus(x.*trcx, rho.*kx, h, [1,2], BCD);

% Boundary conditions
bnd_TRC = zeros(size(TRCx));
if ~isnan(trcwall(1)); bnd_TRC = bnd_TRC + ((trcwall(1,:,:).*rho)-TRCx)./(tau_T+dt) .* topshape; end
if ~isnan(trcwall(2)); bnd_TRC = bnd_TRC + ((trcwall(2,:,:).*rho)-TRCx)./(tau_T+dt) .* botshape; end

% Total rate of change
dTRCmdt = - advn_TRCm + diff_TRCm + Gmtrc + Gemt + Gint;         
dTRCxdt = - advn_TRCx + diff_TRCx + Gxtrc + Gext        + bnd_TRC;

dTRCmdt_lim0 = -TRCm/3/dt;
dTRCxdt_lim0 = -TRCx/3/dt;

dTRCmdt = max(dTRCmdt_lim0,dTRCmdt);
dTRCxdt = max(dTRCxdt_lim0,dTRCxdt);

% Time integration
res_TRCm = (a1*TRCm-a2*TRCmo-a3*TRCmoo) - (b1*dTRCmdt + b2*dTRCmdto + b3*dTRCmdtoo)*dt;
res_TRCx = (a1*TRCx-a2*TRCxo-a3*TRCxoo) - (b1*dTRCxdt + b2*dTRCxdto + b3*dTRCxdtoo)*dt;

[trcm,GHST.TRCm,FHST.TRCm,specrad.TRCm] = iterate(trcm,res_TRCm./max(mutra.*rho,M),specrad.TRCm,GHST.TRCm,FHST.TRCm,itpar,iter);
[trcx,GHST.TRCx,FHST.TRCx,specrad.TRCx] = iterate(trcx,res_TRCx./max(mutra.*rho,X),specrad.TRCx,GHST.TRCx,FHST.TRCx,itpar,iter);

% impose min/max limits on component densities
trcm = max(0,trcm);
trcx = max(0,trcx);

trc  = m.*trcm + x.*trcx;
TRCm = M.*trcm;
TRCx = X.*trcx;
TRC  = TRCm + TRCx;
