%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

tic;

%***  update heat content (entropy) density

% heat advection
[advn_Sm,qz_advn_Sm,qx_advn_Sm] = advect(M.*sm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);  % melt  advection
[advn_Sx,qz_advn_Sx,qx_advn_Sx] = advect(X.*sx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);  % solid advection

[diff_S ,qz_diff_S ,qx_diff_S ] = diffus(T ,        kT./T,h,[1,2],BCD);
[diff_Sd,qz_diff_Sd,qx_diff_Sd] = diffus(Tp,M.*cPm.*kd./T,h,[1,2],BCD);
diff_S    = diff_S + diff_Sd;
qz_diff_S = qz_diff_S + qz_diff_Sd;
qx_diff_S = qx_diff_S + qx_diff_Sd;

% heat dissipation
diss_h = diss ./ T;

% boundary layers
bnd_T = zeros(size(S));
if ~isnan(Twall(1)); bnd_T = bnd_T + ((Twall(1,:)+273.15)-T)./(tau_T+dt) .* topshape; end
if ~isnan(Twall(2)); bnd_T = bnd_T + ((Twall(2,:)+273.15)-T)./(tau_T+dt) .* botshape; end
bnd_S = rho.*cP.*bnd_T ./ T;

% total rate of change
dSdt  = - advn_Sm - advn_Sx + diff_S + diss_h + bnd_S + Gems + Gexs + Gins;

% residual of entropy evolution
res_S = (a1*S-a2*So-a3*Soo) - (b1*dSdt + b2*dSdto + b3*dSdtoo)*dt;

% semi-implicit update of bulk entropy density
[S,GHST.S,FHST.S,specrad.S] = iterate(S,res_S,specrad.S,GHST.S,FHST.S,itpar,iter);

% convert entropy S to natural temperature T and potential temperature Tp
s = S./RHO;
[Tp,~ ] = StoT(Tp,s,cat(3,Pt,Pt)*0+Pref,cat(3,m,x),[cPm;cPx],[aTm;aTx],[bPm;bPx],cat(3,rhom0,rhox0),[sref+Dsm;sref],Tref,Pref);
[T ,si] = StoT(T ,s,cat(3,Pt,Pt)       ,cat(3,m,x),[cPm;cPx],[aTm;aTx],[bPm;bPx],cat(3,rhom0,rhox0),[sref+Dsm;sref],Tref,Pref);
sm = si(:,:,1); sx = si(:,:,2);  % read out phase entropies


%***  update major component densities

% major component advection
[advn_Cm,qz_advn_Cm,qx_advn_Cm] = advect(M.*cm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % melt  advection
[advn_Cx,qz_advn_Cx,qx_advn_Cx] = advect(X.*cx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % solid advection

% major component dispersion
[diff_Cm,qz_diff_Cm,qx_diff_Cm] = diffus(cm,M.*kd,h,[1,2],BCD); 
[diff_Cx,qz_diff_Cx,qx_diff_Cx] = diffus(cx,X.*kx,h,[1,2],BCD);

bnd_C = zeros(size(C));
if ~isnan(cwall(1)); bnd_C = bnd_C + ((cwall(1,:,:).*rho)-C)./(tau_T+dt) .* topshape; end
if ~isnan(cwall(2)); bnd_C = bnd_C + ((cwall(2,:,:).*rho)-C)./(tau_T+dt) .* botshape; end

% total rate of change
dCdt = - advn_Cm - advn_Cx + diff_Cm + diff_Cx + bnd_C + Gemc + Gexc + Ginc + bnd_C;                         
  
% residual of major component evolution
res_C = (a1*C-a2*Co-a3*Coo) - (b1*dCdt + b2*dCdto + b3*dCdtoo)*dt;

% semi-implicit update of major component density
[C,GHST.C,FHST.C,specrad.C] = iterate(C,res_C,specrad.C,GHST.C,FHST.C,itpar,iter);

% impose min/max limits on component densities
C = max(0,min(rho, C ));

% convert component density to concentration
RHOC = sum(C,3);
c = C./RHOC;


%*** update phase equilibrium if full reactive coupling on

if Rcouple || frst; phseql; end


%***  update phase fraction densities

% phase advection rates
[advn_X,qz_advn_X,qx_advn_X] = advect(X,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);
[advn_M,qz_advn_M,qx_advn_M] = advect(M,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);

% total rates of change
dXdt   = - advn_X + Gx + Gex;
dMdt   = - advn_M + Gm + Gem + Gin;
drhodt = dXdt + dMdt;

% residual of phase density evolution
res_X = (a1*X-a2*Xo-a3*Xoo) - (b1*dXdt + b2*dXdto + b3*dXdtoo)*dt;
res_M = (a1*M-a2*Mo-a3*Moo) - (b1*dMdt + b2*dMdto + b3*dMdtoo)*dt;

% semi-implicit update of phase fraction densities
res_PHS = cat(3,res_X,res_M);
PHS     = cat(3,X,M);
[PHS,GHST.PHS,FHST.PHS,specrad.PHS] = iterate(PHS,res_PHS,specrad.PHS,GHST.PHS,FHST.PHS,itpar,iter);

% impose min/max limits on phase densities
X     = max(0,min(rho, PHS(:,:,1) ));
M     = max(0,min(rho, PHS(:,:,2) ));


%***  update phase fractions and component concentrations

% update phase fractions
RHO = X + M;
x   = X./RHO; 
m   = M./RHO;

hasx = x >= eps^0.5;
hasm = m >= eps^0.5;


% update major component phase composition
Kx      = reshape(cal.Kx,Nz,Nx,cal.ncmp);
subsol  = m<=eps^0.5 & T<=reshape(cal.Tsol+273.15,Nz,Nx);
supliq  = x<=eps^0.5 & T>=reshape(cal.Tliq+273.15,Nz,Nx);
subsolc = repmat(subsol,1,1,cal.ncmp);
supliqc = repmat(supliq,1,1,cal.ncmp);
rnorm   = 1;  tol  = atol;
it      = 1;  mxit = 100;
upd_cm  = 0.*cm;  upd_cx = 0.*cx;
cm = cmq;  cx = cxq;
while rnorm>tol && it<mxit

    Kx = cx./(cm+eps);

    cm  = c    ./(m + x.*Kx + eps); 
    cx  = c.*Kx./(m + x.*Kx + eps); 

    cm  = cm./sum(cm,3);
    cx  = cx./sum(cx,3);

    r = x.*cx + m.*cm - c;
    r(subsolc) = 0; r(supliqc) = 0;
    rnorm = norm(r(:))./norm(c(:));

    it  = it+1;
end

if (it==mxit && rnorm>tol)
    disp(['!!! Lever rule adjustment not converged after ',num2str(mxit),' iterations !!!']);
end

% fix subsolidus and superliquidus conditions
cx(subsolc) = c(subsolc); x(subsol) = 1; m(subsol) = 0;
cm(supliqc) = c(supliqc); m(supliq) = 1; x(supliq) = 0;

% record timing
TCtime = TCtime + toc - eqtime.*Rcouple;
