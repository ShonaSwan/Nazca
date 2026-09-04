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
[advn_Cm,qz_advn_Cm,qx_advn_Cm] = advect(Cm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % melt  advection
[advn_Cx,qz_advn_Cx,qx_advn_Cx] = advect(Cx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % solid advection

% major component dispersion
[diff_Cm,qz_diff_Cm,qx_diff_Cm] = diffus(Cm./RHO, RHO.*kd, h, [1,2], BCD); 
[diff_Cx,qz_diff_Cx,qx_diff_Cx] = diffus(Cx./RHO, RHO.*kd, h, [1,2], BCD);

bnd_C = zeros(size(C));
if ~isnan(cwall(1)); bnd_C = bnd_C + ((cwall(1,:,:).*RHO)-Cx)./(tau_T+dt) .* topshape; end
if ~isnan(cwall(2)); bnd_C = bnd_C + ((cwall(2,:,:).*RHO)-Cx)./(tau_T+dt) .* botshape; end

% total rate of change
dCmdt = -advn_Cm + diff_Cm + Gmc + Gemc + Ginc;
dCxdt = -advn_Cx + diff_Cx + Gxc + Gexc        + bnd_C;

% residual of major component evolution
res_Cm = (a1*Cm-a2*Cmo-a3*Cmoo) - (b1*dCmdt + b2*dCmdto + b3*dCmdtoo)*dt;
res_Cx = (a1*Cx-a2*Cxo-a3*Cxoo) - (b1*dCxdt + b2*dCxdto + b3*dCxdtoo)*dt;

% semi-implicit update of major component density
[Cm, GHST.Cm, FHST.Cm, specrad.Cm] = iterate(Cm, res_Cm, specrad.Cm, GHST.Cm, FHST.Cm, itpar,iter, [0,nan]);
[Cx, GHST.Cx, FHST.Cx, specrad.Cx] = iterate(Cx, res_Cx, specrad.Cx, GHST.Cx, FHST.Cx, itpar,iter, [0,nan]);

% impose min/max limits on component densities
Cm = max(0,Cm);
Cx = max(0,Cx);

% update phase fractions
M   = sum(Cm,3);
X   = sum(Cx,3);
RHO = M + X;

dMdt = sum(dCmdt,3);
dXdt = sum(dCxdt,3);
drhodt = dMdt + dXdt;

m   = M./RHO;
x   = X./RHO;

C  = Cm + Cx;
c  = C./RHO;

% blend diseq/eq phase and mixture major components
Kx = (cxq+eps)./(cmq+eps);
for i=1:20
    cm  = (1-f_diseq).*c./(m + x.*Kx) + f_diseq.*Cm./max(M,1e-12.*RHO);  cm = cm./sum(cm+eps,3);
    cx  = (1-f_diseq).*c./(m./Kx + x) + f_diseq.*Cx./max(X,1e-12.*RHO);  cx = cx./sum(cx+eps,3);
    Kx  = (cx+eps)./(cm+eps);
    res = c - (m.*cm + x.*cx);
end
Cm = cm.*M;
Cx = cx.*X;

% unpack major and trace elements
majm = cm(:,:,1:cal.ncmp);  majm = majm./sum(majm,3);
majx = cx(:,:,1:cal.ncmp);  majx = majx./sum(majx,3);
maj  = c (:,:,1:cal.ncmp);  maj  = maj ./sum(maj ,3);
trcm = cm(:,:,cal.ncmp+1:end)/trcscale;
trcx = cx(:,:,cal.ncmp+1:end)/trcscale;
trc  = c (:,:,cal.ncmp+1:end)/trcscale;

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
