%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

tic;

%***  update heat content (entropy) density

% heat advection
[advn_Sm,qz_advn_Sm,qx_advn_Sm] = advect(Sm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);  % melt  advection
[advn_Sx,qz_advn_Sx,qx_advn_Sx] = advect(Sx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);  % solid advection

% major component dispersion
[diff_S ,qz_diff_S ,qx_diff_S ] = diffus(T ,        kT./T,h,[1,2],BCD);
[diff_Sm,qz_diff_Sm,qx_diff_Sm] = diffus(Sm./RHO, RHO.*kd, h, [1,2], BCD); 
[diff_Sx,qz_diff_Sx,qx_diff_Sx] = diffus(Sx./RHO, RHO.*kd, h, [1,2], BCD);

diff_S    = diff_S + diff_Sm + diff_Sx;
qz_diff_S = qz_diff_S + qz_diff_Sm + qz_diff_Sx;
qx_diff_S = qx_diff_S + qx_diff_Sm + qx_diff_Sx;

% heat dissipation
diss_h = diss ./ T;

% boundary layers
bnd_T = zeros(size(S));
if ~isnan(Twall(1)); bnd_T = bnd_T + ((Twall(1,:)+273.15)-T)./(tau_T+dt) .* topshape; end
if ~isnan(Twall(2)); bnd_T = bnd_T + ((Twall(2,:)+273.15)-T)./(tau_T+dt) .* botshape; end
bnd_S = RHO.*cP.*bnd_T ./ T;

% total rate of change
dSdt  = - advn_Sm - advn_Sx + diff_S + diss_h + Gems + Gexs + Gins + bnd_S;

% residual of entropy evolution
res_S = (a1*S-a2*So-a3*Soo) - (b1*dSdt + b2*dSdto + b3*dSdtoo)*dt;

% semi-implicit update of bulk entropy density
[S,GHST.S,FHST.S,specrad.S] = iterate(S,res_S,specrad.S,GHST.S,FHST.S,itpar,iter,[nan,nan]);


%***  update major and trace components

% major component advection
[advn_Cm,qz_advn_Cm,qx_advn_Cm] = advect(CCm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % melt  advection
[advn_Cx,qz_advn_Cx,qx_advn_Cx] = advect(CCx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % solid advection

% major component dispersion
[diff_Cm,qz_diff_Cm,qx_diff_Cm] = diffus(CCm./RHO, RHO.*kd, h, [1,2], BCD); 
[diff_Cx,qz_diff_Cx,qx_diff_Cx] = diffus(CCx./RHO, RHO.*kd, h, [1,2], BCD);

bnd_C = zeros(size(CC));
if ~isnan(cwall(1)); bnd_C = bnd_C + ((ccwall(1,:,:).*RHO)-CCx)./(tau_T+dt) .* topshape; end
if ~isnan(cwall(2)); bnd_C = bnd_C + ((ccwall(2,:,:).*RHO)-CCx)./(tau_T+dt) .* botshape; end

% total rate of change
dCmdt = -advn_Cm + diff_Cm + Gmc + Gemc + Ginc;
dCxdt = -advn_Cx + diff_Cx + Gxc + Gexc        + bnd_C;

% residual of major component evolution
res_Cm = (a1*CCm-a2*CCmo-a3*CCmoo) - (b1*dCmdt + b2*dCmdto + b3*dCmdtoo)*dt;
res_Cx = (a1*CCx-a2*CCxo-a3*CCxoo) - (b1*dCxdt + b2*dCxdto + b3*dCxdtoo)*dt;

% semi-implicit update of major component density
[CCm, GHST.Cm, FHST.Cm, specrad.Cm] = iterate(CCm, res_Cm, specrad.Cm, GHST.Cm, FHST.Cm, itpar,iter, [0,nan]);
[CCx, GHST.Cx, FHST.Cx, specrad.Cx] = iterate(CCx, res_Cx, specrad.Cx, GHST.Cx, FHST.Cx, itpar,iter, [0,nan]);

% impose min/max limits on component densities
CCm = max(0,CCm);
CCx = max(0,CCx);

% update phase fractions
M   = sum(CCm,3);
X   = sum(CCx,3);
RHO = M + X;

dMdt = sum(dCmdt,3);
dXdt = sum(dCxdt,3);
drhodt = dMdt + dXdt;

m   = M./RHO;
x   = X./RHO;

CC  = CCm + CCx;
cc  = CC./RHO;

% blend diseq/eq phase and mixture major components
Kx = (ccxq+eps)./(ccmq+eps);
for i=1:20
    ccm = (1-f_diseq).*cc./(m + x.*Kx) + f_diseq.*CCm./max(M,1e-9.*RHO);  ccm = ccm./sum(ccm+eps,3);
    ccx = (1-f_diseq).*cc./(m./Kx + x) + f_diseq.*CCx./max(X,1e-9.*RHO);  ccx = ccx./sum(ccx+eps,3);
    Kx  = (ccx+eps)./(ccm+eps);
    res = cc - (m.*ccm + x.*ccx);
end
CCm = ccm.*M;
CCx = ccx.*X;
CC  = CCm + CCx;

% unpack major and trace elements
cm = ccm(:,:,1:cal.ncmp);  cm = cm./sum(cm,3);
cx = ccx(:,:,1:cal.ncmp);  cx = cx./sum(cx,3);
c  = cc (:,:,1:cal.ncmp);  c  = c ./sum(c ,3);
trcm = ccm(:,:,cal.ncmp+1:end)/trcscale;
trcx = ccx(:,:,cal.ncmp+1:end)/trcscale;
trc  = cc (:,:,cal.ncmp+1:end)/trcscale;

% truncate residuals to diseq/eq domain
res_Cm = f_diseq.*res_Cm;
res_Cx = f_diseq.*res_Cx;

% convert entropy S to natural temperature T and potential temperature Tp
s = S./RHO;
[Tp,~ ] = StoT(Tp,s,cat(3,Pt,Pt)*0+Pref,cat(3,m,x),[cPm;cPx],[aTm;aTx],[bPm;bPx],cat(3,rhom0,rhox0),[sref+Dsm;sref],Tref,Pref);
[T ,si] = StoT(T ,s,cat(3,Pt,Pt)       ,cat(3,m,x),[cPm;cPx],[aTm;aTx],[bPm;bPx],cat(3,rhom0,rhox0),[sref+Dsm;sref],Tref,Pref);
sm = si(:,:,1); sx = si(:,:,2);  % read out phase entropies
Sm = sm.*M;
Sx = sx.*X;


%*** update phase equilibrium if full reactive coupling on

if Rcouple; phseql; end

% record timing
TCtime = TCtime + toc - eqtime.*Rcouple;
