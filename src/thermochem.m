%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

tic;

%***  update heat content (entropy) density

% heat advection
[advn_Sm,qz_advn_Sm,qx_advn_Sm] = advect(Sm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);  % melt  advection
[advn_Sx,qz_advn_Sx,qx_advn_Sx] = advect(Sx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);  % solid advection

% major component dispersion
[diff_S ,qz_diff_S ,qx_diff_S ] = diffus(T ,        kT./T,h,[1,2],BCD);
[diff_Sm,qz_diff_Sm,qx_diff_Sm] = diffus(m.*sm, rho.*kd, h, [1,2], BCD); 
[diff_Sx,qz_diff_Sx,qx_diff_Sx] = diffus(x.*sx, rho.*kd, h, [1,2], BCD);

diff_S    = diff_S + diff_Sm + diff_Sx;
qz_diff_S = qz_diff_S + qz_diff_Sm + qz_diff_Sx;
qx_diff_S = qx_diff_S + qx_diff_Sm + qx_diff_Sx;

% heat dissipation
diss_h = diss ./ T;

% boundary layers
bnd_T = zeros(size(S));
if ~isnan(Twall(1)); bnd_T = bnd_T + ((Twall(1,:)+273.15)-T)./(tau_T+dt) .* topshape; end
if ~isnan(Twall(2)); bnd_T = bnd_T + ((Twall(2,:)+273.15)-T)./(tau_T+dt) .* botshape; end
bnd_S = rho.*cP.*bnd_T ./ T;

% total rate of change
dSdt  = - advn_Sm - advn_Sx + diff_S + diss_h + Gems + Gexs + Gins + bnd_S;

% residual of entropy evolution
res_S = (a1*S-a2*So-a3*Soo) - (b1*dSdt + b2*dSdto + b3*dSdtoo)*dt;

% semi-implicit update of bulk entropy density
[S,GHST.S,FHST.S,specrad.S] = iterate(S,res_S,specrad.S,GHST.S,FHST.S,itpar,iter,[nan,nan]);


%***  update major components

% major component advection
[advn_Cm,qz_advn_Cm,qx_advn_Cm] = advect(Cm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % melt  advection
[advn_Cx,qz_advn_Cx,qx_advn_Cx] = advect(Cx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % solid advection

% major component dispersion
[diff_Cm,qz_diff_Cm,qx_diff_Cm] = diffus(m.*cm, rho.*kd, h, [1,2], BCD); 
[diff_Cx,qz_diff_Cx,qx_diff_Cx] = diffus(x.*cx, rho.*kd, h, [1,2], BCD);

bnd_C = zeros(size(C));
if ~isnan(cwall(1)); bnd_C = bnd_C + ((cwall(1,:,:).*rho)-C)./(tau_T+dt) .* topshape; end
if ~isnan(cwall(2)); bnd_C = bnd_C + ((cwall(2,:,:).*rho)-C)./(tau_T+dt) .* botshape; end

% total rate of change
dCmdt = -advn_Cm + diff_Cm + Gmc + Gemc + Ginc;
dCxdt = -advn_Cx + diff_Cx + Gxc + Gexc        + bnd_C;
dCdt  = dCmdt + dCxdt;

% residual of major component evolution
res_Cm = (a1*Cm-a2*Cmo-a3*Cmoo) - (b1*dCmdt + b2*dCmdto + b3*dCmdtoo)*dt;
res_Cx = (a1*Cx-a2*Cxo-a3*Cxoo) - (b1*dCxdt + b2*dCxdto + b3*dCxdtoo)*dt;
res_C  = (a1*C -a2*Co -a3*Coo ) - (b1*dCdt  + b2*dCdto  + b3*dCdtoo )*dt;

% semi-implicit update of major component density
[Cm, GHST.Cm, FHST.Cm, specrad.Cm] = iterate(Cm, res_Cm, specrad.Cm, GHST.Cm, FHST.Cm, itpar,iter, [0,nan]);
[Cx, GHST.Cx, FHST.Cx, specrad.Cx] = iterate(Cx, res_Cx, specrad.Cx, GHST.Cx, FHST.Cx, itpar,iter, [0,nan]);
[C , GHST.C , FHST.C , specrad.C ] = iterate(C , res_C , specrad.C , GHST.C , FHST.C , itpar,iter, [0,nan]);

% impose min/max limits on component densities
Cm = max(0,min(rho,Cm));  %cm = cm./sum(cm+eps,3);
Cx = max(0,min(rho,Cx));  %cx = cx./sum(cx+eps,3);
C  = max(0,min(rho,C ));  %c  = c ./sum(c +eps,3);


%***  update trace components

% trace component advection
[advn_TRCm,qz_advn_TRCm,qx_advn_TRCm] = advect(TRCm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % melt  advection
[advn_TRCx,qz_advn_TRCx,qx_advn_TRCx] = advect(TRCx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % solid advection

% trace component dispersion
[diff_TRCm,qz_diff_TRCm,qx_diff_TRCm] = diffus(m.*trcm, rho.*kd, h, [1,2], BCD); 
[diff_TRCx,qz_diff_TRCx,qx_diff_TRCx] = diffus(x.*trcx, rho.*kx, h, [1,2], BCD);

% boundary conditions
bnd_TRC = zeros(size(TRCx));
if ~isnan(trcwall(1)); bnd_TRC = bnd_TRC + ((trcwall(1,:,:).*rho)-TRCx)./(tau_T+dt) .* topshape; end
if ~isnan(trcwall(2)); bnd_TRC = bnd_TRC + ((trcwall(2,:,:).*rho)-TRCx)./(tau_T+dt) .* botshape; end

% total rate of change
dTRCmdt = - advn_TRCm + diff_TRCm + Gmtrc + Gemt + Gint;         
dTRCxdt = - advn_TRCx + diff_TRCx + Gxtrc + Gext        + bnd_TRC;
dTRCdt  = dTRCmdt + dTRCxdt;

% time integration
res_TRCm = (a1*TRCm-a2*TRCmo-a3*TRCmoo) - (b1*dTRCmdt + b2*dTRCmdto + b3*dTRCmdtoo)*dt;
res_TRCx = (a1*TRCx-a2*TRCxo-a3*TRCxoo) - (b1*dTRCxdt + b2*dTRCxdto + b3*dTRCxdtoo)*dt;
res_TRC  = (a1*TRC -a2*TRCo -a3*TRCoo ) - (b1*dTRCdt  + b2*dTRCdto  + b3*dTRCdtoo )*dt;

[TRCm,GHST.TRCm,FHST.TRCm,specrad.TRCm] = iterate(TRCm,res_TRCm,specrad.TRCm,GHST.TRCm,FHST.TRCm,itpar,iter,[0,nan]);
[TRCx,GHST.TRCx,FHST.TRCx,specrad.TRCx] = iterate(TRCx,res_TRCx,specrad.TRCx,GHST.TRCx,FHST.TRCx,itpar,iter,[0,nan]);
[TRC ,GHST.TRC ,FHST.TRC ,specrad.TRC ] = iterate(TRC ,res_TRC ,specrad.TRC ,GHST.TRC ,FHST.TRC ,itpar,iter,[0,nan]);

% impose min/max limits on component densities
TRCm = max(0,TRCm);
TRCx = max(0,TRCx);


%***  update phase fractions

% phase advection rates
[advn_X,qz_advn_X,qx_advn_X] = advect(X,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);
[advn_M,qz_advn_M,qx_advn_M] = advect(M,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);

% phase dispersion
[diff_M,qz_diff_M,qx_diff_M] = diffus(m, rho.*kd, h, [1,2], BCD); 
[diff_X,qz_diff_X,qx_diff_X] = diffus(x, rho.*kd, h, [1,2], BCD);

% total rates of change
dMdt = - advn_M + diff_M + Gm + Gem + Gin;
dXdt = - advn_X + diff_X + Gx + Gex      ;
drhodt = - advn_M - advn_X + Gem + Gin + Gex;

% residual of phase density evolution
res_X = (a1*X-a2*Xo-a3*Xoo) - (b1*dXdt + b2*dXdto + b3*dXdtoo)*dt;
res_M = (a1*M-a2*Mo-a3*Moo) - (b1*dMdt + b2*dMdto + b3*dMdtoo)*dt;

[X,GHST.X,FHST.X,specrad.X] = iterate(X,res_X,specrad.X,GHST.X,FHST.X,itpar,iter, [0,nan]);
[M,GHST.M,FHST.M,specrad.M] = iterate(M,res_M,specrad.M,GHST.M,FHST.M,itpar,iter, [0,nan]);

% impose min/max limits and normalise component densities
M  = max(0,min(rho,M));
X  = max(0,min(rho,X));

m  = M./(M+X);
x  = X./(M+X);


%***  update partial densities

% partial phase densities
% X  = x.*rho;
% M  = m.*rho;

% blend diseq/eq phase and mixture major components
cq  = C./(M+X);  cq = cq./sum(cq+eps,3);
cmd = (1-f_diseq).*cmq + f_diseq .* Cm./max(M,1e-9.*rho);  cmd = cmd./sum(cmd+eps,3);
cxd = (1-f_diseq).*cxq + f_diseq .* Cx./max(X,1e-9.*rho);  cxd = cxd./sum(cxd+eps,3);
c   = (1-f_diseq).*cq  + f_diseq .* (m.*cmd + x.*cxd);
Kx  = (cxq+eps)./(cmq+eps);
for i=1:20
cm = (1-f_diseq).*c./(m + x.*Kx) + f_diseq.*cmd;  cm = cm./sum(cm+eps,3);
cx = (1-f_diseq).*c./(m./Kx + x) + f_diseq.*cxd;  cx = cx./sum(cx+eps,3);
Kx = (cx+eps)./(cm+eps);
end

% partial major component densities
% Cm = cm.*M;
% Cx = cx.*X;
% C  = Cm + Cx;

% truncate residuals to diseq/eq domain
res_Cm = f_diseq.*res_Cm;
res_Cx = f_diseq.*res_Cx;
res_C  = (1-f_diseq).*res_C;

% blend diseq/eq phase and mixture major components
trcq  = TRC./(M+X);
trcmd = (1-f_diseq).*trcmq + f_diseq .* TRCm./max(M,1e-9.*rho);
trcxd = (1-f_diseq).*trcxq + f_diseq .* TRCx./max(X,1e-9.*rho);
trc   = (1-f_diseq).*trcq  + f_diseq .* (m.*trcmd + x.*trcxd);
trcm  = (1-f_diseq).*trc./(m + x.*Ktrc) + f_diseq.*trcm;
trcx  = (1-f_diseq).*trc./(m./Ktrc + x) + f_diseq.*trcx;

% partial trace component densities
% TRCm = M.*trcm;
% TRCx = X.*trcx;
% TRC  = TRCm + TRCx;

% convert entropy S to natural temperature T and potential temperature Tp
s = S./(M+X);
[Tp,~ ] = StoT(Tp,s,cat(3,Pt,Pt)*0+Pref,cat(3,m,x),[cPm;cPx],[aTm;aTx],[bPm;bPx],cat(3,rhom0,rhox0),[sref+Dsm;sref],Tref,Pref);
[T ,si] = StoT(T ,s,cat(3,Pt,Pt)       ,cat(3,m,x),[cPm;cPx],[aTm;aTx],[bPm;bPx],cat(3,rhom0,rhox0),[sref+Dsm;sref],Tref,Pref);
sm = si(:,:,1); sx = si(:,:,2);  % read out phase entropies

% partial entropy densities
% Sm = M.*sm;
% Sx = X.*sx;
% S  = Sm + Sx;


%*** update phase equilibrium if full reactive coupling on

if Rcouple; phseql; end

% record timing
TCtime = TCtime + toc - eqtime.*Rcouple;
