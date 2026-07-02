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
dTRCmdt = - advn_TRCm + diff_TRCm + Gemt + Gint + Gmtrc + Gmtrcex;              % re-equilibration
dTRCxdt = - advn_TRCx + diff_TRCx + Gext        + Gxtrc + Gxtrcex + bnd_TRC;    % re-equilibration

% Time integration
res_TRCm = (a1*TRCm-a2*TRCmo-a3*TRCmoo) - (b1*dTRCmdt + b2*dTRCmdto + b3*dTRCmdtoo)*dt;
res_TRCx = (a1*TRCx-a2*TRCxo-a3*TRCxoo) - (b1*dTRCxdt + b2*dTRCxdto + b3*dTRCxdtoo)*dt;

[TRCm,GHST.TRCm,FHST.TRCm,specrad.TRCm] = iterate(TRCm,res_TRCm,specrad.TRCm,GHST.TRCm,FHST.TRCm,itpar,iter);
[TRCx,GHST.TRCx,FHST.TRCx,specrad.TRCx] = iterate(TRCx,res_TRCx,specrad.TRCx,GHST.TRCx,FHST.TRCx,itpar,iter);

trcm = TRCm ./ max(1e-3,M); 
trcx = TRCx ./ max(1e-3,X);

% TRC = TRCm + TRCx;

