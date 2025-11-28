%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

tic;

%***  update heat content (entropy) density

% heat advection
[advn_Sm,qz_advn_Sm,qx_advn_Sm] = advect(M.*sm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);  % melt  advection
[advn_Sx,qz_advn_Sx,qx_advn_Sx] = advect(X.*sx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);  % solid advection

[dffn_S,qz_dffn_S,qx_dffn_S] = diffus(T,kT./T,h,[1,2],BCD);

% heat dissipation
diss_h = diss ./ T;

% boundary layers
bnd_T = zeros(size(S));
if ~isnan(Twall(1)); bnd_T = bnd_T + ((Twall(1)+273.15)-T)./(tau_T+dt) .* topshape; end
if ~isnan(Twall(2)); bnd_T = bnd_T + ((Twall(2)+273.15)-T)./(tau_T+dt) .* botshape; end
bnd_S = RHO.*cP.*bnd_T ./ T;


% total rate of change
dSdt  = - advn_Sm - advn_Sx + dffn_S + diss_h + bnd_S + Gems + Gexs + Gins;

% residual of entropy evolution
res_S = (a1*S-a2*So-a3*Soo)/dt - (b1*dSdt + b2*dSdto + b3*dSdtoo);

% semi-implicit update of bulk entropy density
[S,GHST.S,FHST.S,specrad.S] = iterate(S,res_S*dt/a1,specrad.S,GHST.S,FHST.S,itpar,iter*~frst);

% convert entropy S to natural temperature T and potential temperature Tp
[Tp,~ ] = StoT(Tp,S./RHO,cat(3,Pt,Ptx)*0+Pref,cat(3,m,x),[cPm;cPx],[aTm;aTx],[bPm;bPx],cat(3,rhom0,rhox0),[sref+Dsm;sref],Tref,Pref);
[T ,si] = StoT(T ,S./RHO,cat(3,Pt,Ptx)       ,cat(3,m,x),[cPm;cPx],[aTm;aTx],[bPm;bPx],cat(3,rhom0,rhox0),[sref+Dsm;sref],Tref,Pref);
sm = si(:,:,1); sx = si(:,:,2);  % read out phase entropies


%***  update major component densities

% major component advection
[advn_Cm,qz_advn_Cm,qx_advn_Cm] = advect(M.*cm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % melt  advection
[advn_Cx,qz_advn_Cx,qx_advn_Cx] = advect(X.*cx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % solid advection

% total rate of change
dCdt = - advn_Cm - advn_Cx + bnd_C + Gemc + Gexc + Ginc;                                      
  
% residual of major component evolution
res_C = (a1*C-a2*Co-a3*Coo)/dt - (b1*dCdt + b2*dCdto + b3*dCdtoo);

% semi-implicit update of major component density
[C,GHST.C,FHST.C,specrad.C] = iterate(C,res_C*dt/a1,specrad.C,GHST.C,FHST.C,itpar,iter*~frst);

% impose min/max limits on component densities
C = max(0,min(rho, C ));

% convert component density to concentration
c = C./sum(C,3);


%*** update phase equilibrium if full reactive coupling on

if Rcouple; phseql; end


%***  update phase fraction densities

% phase advection rates
[advn_X,qz_advn_X,qx_advn_X] = advect(X,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);
[advn_M,qz_advn_M,qx_advn_M] = advect(M,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);

% total rates of change
dXdt   = - advn_X + Gx + Gex;
dMdt   = - advn_M + Gm + Gem + Gin;
drhodt = - advn_X - advn_M + Gem + Gex + Gin;

% residual of phase density evolution
res_X = (a1*X-a2*Xo-a3*Xoo)/dt - (b1*dXdt + b2*dXdto + b3*dXdtoo);
res_M = (a1*M-a2*Mo-a3*Moo)/dt - (b1*dMdt + b2*dMdto + b3*dMdtoo);

% semi-implicit update of phase fraction densities
res_PHS = cat(3,res_X,res_M);
PHS     = cat(3,X,M);
[PHS,GHST.PHS,FHST.PHS,specrad.PHS] = iterate(PHS,res_PHS*dt/a1,specrad.PHS,GHST.PHS,FHST.PHS,itpar,iter*~frst);

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
cx(subsolc) = cxq(subsolc); x(subsol) = xq(subsol); m(subsol) = 0;
cm(supliqc) = cmq(supliqc); m(supliq) = mq(supliq); x(supliq) = 0;

% record timing
TCtime = TCtime + toc - eqtime.*Rcouple;
