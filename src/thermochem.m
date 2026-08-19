%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

tic;

%***  update heat content (entropy) density

% heat advection
[advn_Sm,qz_advn_Sm,qx_advn_Sm] = advect(M.*sm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);  % melt  advection
[advn_Sx,qz_advn_Sx,qx_advn_Sx] = advect(X.*sx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);  % solid advection

% major component dispersion
[diff_S ,qz_diff_S ,qx_diff_S ] = diffus(T ,        kT./T,h,[1,2],BCD);
[diff_Sm,qz_diff_Sm,qx_diff_Sm] = diffus(m.*sm, rho.*kd, h, [1,2], BCD); 
[diff_Sx,qz_diff_Sx,qx_diff_Sx] = diffus(x.*sx, rho.*kx, h, [1,2], BCD);

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
[S,GHST.S,FHST.S,specrad.S] = iterate(S,res_S,specrad.S,GHST.S,FHST.S,itpar,iter);
% S(1,:)     = (S(1,:) + cPx.*rhox(1,:).*log((T0+273.15)./Tref))/2;
% ref_S(1,:) = 0;

% convert entropy S to natural temperature T and potential temperature Tp
s = S./RHO;
[Tp,~ ] = StoT(Tp,s,cat(3,Pt,Pt)*0+Pref,cat(3,m,x),[cPm;cPx],[aTm;aTx],[bPm;bPx],cat(3,rhom0,rhox0),[sref+Dsm;sref],Tref,Pref);
[T ,si] = StoT(T ,s,cat(3,Pt,Pt)       ,cat(3,m,x),[cPm;cPx],[aTm;aTx],[bPm;bPx],cat(3,rhom0,rhox0),[sref+Dsm;sref],Tref,Pref);
sm = si(:,:,1); sx = si(:,:,2);  % read out phase entropies

%% 
%***  update major component densities

% major component advection
[advn_Cm,qz_advn_Cm,qx_advn_Cm] = advect(Cm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % melt  advection
[advn_Cx,qz_advn_Cx,qx_advn_Cx] = advect(Cx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % solid advection

% major component dispersion
[diff_Cm,qz_diff_Cm,qx_diff_Cm] = diffus(m.*cm, rho.*kd, h, [1,2], BCD); 
[diff_Cx,qz_diff_Cx,qx_diff_Cx] = diffus(x.*cx, rho.*kx, h, [1,2], BCD);

bnd_C = zeros(size(C));
if ~isnan(cwall(1)); bnd_C = bnd_C + ((cwall(1,:,:).*rho)-C)./(tau_T+dt) .* topshape; end
if ~isnan(cwall(2)); bnd_C = bnd_C + ((cwall(2,:,:).*rho)-C)./(tau_T+dt) .* botshape; end

% total rate of change
dCmdt = -advn_Cm + diff_Cm + Gmc + Gemc + Ginc;
dCxdt = -advn_Cx + diff_Cx + Gxc + Gexc        + bnd_C;

dCmdt_lim0 = -Cm/3/dt;
dCxdt_lim0 = -Cx/3/dt;

dCmdt_lim1 =  (rho-Cm)/3/dt;
dCxdt_lim1 =  (rho-Cx)/3/dt;

dCmdt = max(dCmdt_lim0,min(dCmdt_lim1,dCmdt));
dCxdt = max(dCxdt_lim0,min(dCxdt_lim1,dCxdt));

% residual of major component evolution
res_Cm = (a1*Cm-a2*Cmo-a3*Cmoo) - (b1*dCmdt + b2*dCmdto + b3*dCmdtoo)*dt;
res_Cx = (a1*Cx-a2*Cxo-a3*Cxoo) - (b1*dCxdt + b2*dCxdto + b3*dCxdtoo)*dt;

% semi-implicit update of major component density
[cm, GHST.Cm, FHST.Cm, specrad.Cm] = iterate(cm, res_Cm./max(mutra.*rho,M), specrad.Cm, GHST.Cm, FHST.Cm, itpar, iter);
[cx, GHST.Cx, FHST.Cx, specrad.Cx] = iterate(cx, res_Cx./max(mutra.*rho,X), specrad.Cx, GHST.Cx, FHST.Cx, itpar, iter);

% impose min/max limits on component densities
cm = max(0,min(1,cm));  cm = cm./sum(cm,3);
cx = max(0,min(1,cx));  cx = cx./sum(cx,3);

c  = m.*cm + x.*cx;
Cm = cm.*M;
Cx = cx.*X;
C  = Cm + Cx;


%*** update phase equilibrium if full reactive coupling on

if Rcouple || frst; phseql; end


%% 

%***  update phase fraction densities

% phase advection rates
[advn_X,qz_advn_X,qx_advn_X] = advect(X,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);
[advn_M,qz_advn_M,qx_advn_M] = advect(M,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);

% total rates of change
dXdt = - advn_X + Gx + Gem + Gin;
dMdt = - advn_M + Gm + Gex      ;

dMdt_lim0 = - M/3/dt;
dXdt_lim0 = - X/3/dt;

dMdt_lim1 = (rho-M)/3/dt;
dXdt_lim1 = (rho-X)/3/dt;

dMdt = max(dMdt_lim0,min(dMdt_lim1,dMdt));
dXdt = max(dXdt_lim0,min(dXdt_lim1,dXdt));

% dMdt   = sum(dCmdt,3);
% dXdt   = sum(dCxdt,3);

drhodt = dXdt + dMdt;

% residual of phase density evolution
res_X = (a1*X-a2*Xo-a3*Xoo) - (b1*dXdt + b2*dXdto + b3*dXdtoo)*dt;
res_M = (a1*M-a2*Mo-a3*Moo) - (b1*dMdt + b2*dMdto + b3*dMdtoo)*dt;

% semi-implicit update of phase fraction densities
res_PHS = cat(3,res_X,res_M);
phs     = cat(3,x,m);
[phs,GHST.PHS,FHST.PHS,specrad.PHS] = iterate(phs,res_PHS./rho,specrad.PHS,GHST.PHS,FHST.PHS,itpar,iter);

% impose min/max limits on phase densities
phs   = max(0,min(1, phs ));  phs = phs./sum(phs,3);

x = phs(:,:,1);
m = phs(:,:,2);
X = x.*rho;
M = m.*rho;
RHO = X + M;

hasx = x >= eps^0.5;
hasm = m >= eps^0.5;

% Kx      = reshape(cal.Kx,Nz,Nx,cal.ncmp);
% subsol  = m<=eps^0.5 & T<=reshape(cal.Tsol+273.15,Nz,Nx);
% supliq  = x<=eps^0.5 & T>=reshape(cal.Tliq+273.15,Nz,Nx);
% subsolc = repmat(subsol,1,1,cal.ncmp);
% supliqc = repmat(supliq,1,1,cal.ncmp);
% rnorm   = 1;  tol  = atol;
% it      = 1;  mxit = 100;
% upd_cm  = 0.*cm;  upd_cx = 0.*cx;
% cm = cmq;  cx = cxq;
% while rnorm>tol && it<mxit
% 
%     Kx = cx./(cm+eps);
% 
%     cm  = c    ./(m + x.*Kx + eps); 
%     cx  = c.*Kx./(m + x.*Kx + eps); 
% 
%     cm  = cm./sum(cm,3);
%     cx  = cx./sum(cx,3);
% 
%     r = x.*cx + m.*cm - c;
%     r(subsolc) = 0; r(supliqc) = 0;
%     rnorm = norm(r(:))./norm(c(:));
% 
%     it  = it+1;
% end
% 
% if (it==mxit && rnorm>tol)
%     disp(['!!! Lever rule adjustment not converged after ',num2str(mxit),' iterations !!!']);
% end

% fix subsolidus and superliquidus conditions
% cx(subsolc) = c(subsolc); x(subsol) = 1; m(subsol) = 0;
% cm(supliqc) = c(supliqc); m(supliq) = 1; x(supliq) = 0;

% record timing
TCtime = TCtime + toc - eqtime.*Rcouple;
