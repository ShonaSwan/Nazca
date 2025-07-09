%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

tic;

if iter==1; upd_S = 0; upd_C = 0; upd_M = 0; upd_X = 0; end

%***  update heat content (entropy) density

% heat advection
advn_S = - advect(M.*sm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % melt  advection
         - advect(X.*sx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % solid advection

diff_S = diffus(T,kT./T,h,[1,2],BCD);% + diffus(Tp,ks,h,[1,2],BCD);

% heat dissipation (switch this off for now)
diss_h = diss ./ T;

% boundary layers
bnd_T = zeros(size(S));
if ~isnan(Twall(1)); bnd_T = bnd_T + ((Twall(1)+273.15)-T)./(tau_T+dt) .* topshape; end
if ~isnan(Twall(2)); bnd_T = bnd_T + ((Twall(2)+273.15)-T)./(tau_T+dt) .* botshape; end
bnd_S = RHO.*cP.*bnd_T ./ T;


% total rate of change
dSdt  = advn_S + diff_S + diss_h + bnd_S +sm.*Gem + sx.*Gex;

% residual of entropy evolution
res_S = (a1*S-a2*So-a3*Soo)/dt - (b1*dSdt + b2*dSdto + b3*dSdtoo);

% semi-implicit update of bulk entropy density
upd_S = - alpha*res_S*dt/a1 + beta*upd_S;
S     = S + upd_S;

% convert entropy S to natural temperature T and potential temperature Tp
[Tp,~ ] = StoT(Tp,S./RHO,cat(3,Pt,Pt)*0+Pref,cat(3,m,x),[cPm;cPx],[aTm;aTx],[bPm;bPx],cat(3,rhom0,rhox0),[sref;sref+Dsx],Tref,Pref);
[T ,si] = StoT(T ,S./RHO,cat(3,Pt,Pt)       ,cat(3,m,x),[cPm;cPx],[aTm;aTx],[bPm;bPx],cat(3,rhom0,rhox0),[sref;sref+Dsx],Tref,Pref);
sm = si(:,:,1); sx = si(:,:,2);  % read out phase entropies


%***  update major component densities

% major component advection
advn_C = - advect(M.*cm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % melt  advection
         - advect(X.*cx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % solid advection

% major component diffusion (regularisation)
diff_C = diffus(cm,M.*kc,h,[1,2],BCD) + diffus(cx,X.*kc,h,[1,2],BCD);

% boundary layers
bnd_C = zeros(size(C));
for i = 1:cal.ncmp
    if ~isnan(cwall(1)); bnd_C(:,:,i) = bnd_C(:,:,i) + (RHO.*cwall(1,i)-C(:,:,i)).*mu./tau_a .* topshape; end
    if ~isnan(cwall(2)); bnd_C(:,:,i) = bnd_C(:,:,i) + (RHO.*cwall(2,i)-C(:,:,i)).*mu./tau_a .* botshape; end
end

% total rate of change
dCdt = advn_C + diff_C + bnd_C + Gemc + Gexc;                                            
  
% residual of major component evolution
res_C = (a1*C-a2*Co-a3*Coo)/dt - (b1*dCdt + b2*dCdto + b3*dCdtoo);

% semi-implicit update of major component density
upd_C = max(-C, - alpha*res_C*dt/a1 + beta*upd_C );
C     = C + upd_C;

% convert component density to concentration
c = C./sum(C,3);


%*** update phase equilibrium if full reactive coupling on
if Rcouple; phseql; end


%***  update phase fraction densities

% phase advection rates
advn_X   = - advect(X,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_M   = - advect(M,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);
drhodt   = advn_X+advn_M+Gem+Gex;

% total rates of change
dXdt   = advn_X + Gx + Gex;
dMdt   = advn_M + Gm + Gem;

% residual of phase density evolution
res_X = (a1*X-a2*Xo-a3*Xoo)/dt - (b1*dXdt + b2*dXdto + b3*dXdtoo);
res_M = (a1*M-a2*Mo-a3*Moo)/dt - (b1*dMdt + b2*dMdto + b3*dMdtoo);

% semi-implicit update of phase fraction densities
upd_X = max(-X, - alpha*res_X*dt/a1 + beta*upd_X );
upd_M = max(-M, - alpha*res_M*dt/a1 + beta*upd_M );
X     = X + upd_X;
M     = M + upd_M;

%***  update phase fractions and component concentrations

% update phase fractions
RHO = X + M;
x = X./RHO; 
m = M./RHO;

hasx = x >= eps^0.5;
hasm = m >= eps^0.5;

if reactive

% update major component phase composition
Kx      = reshape(cal.Kx,Nz,Nx,cal.ncmp);
%Kf      = reshape(cal.Kf,Nz,Nx,cal.ncmp);
subsol  = m<=eps^0.5 & T<=reshape(cal.Tsol+273.15,Nz,Nx);
supliq  = x<=eps^0.5 & T>=reshape(cal.Tliq+273.15,Nz,Nx);
subsolc = repmat(subsol,1,1,cal.ncmp);
supliqc = repmat(supliq,1,1,cal.ncmp);
rnorm   = 1;  tol  = atol*10;
it      = 1;  mxit = 100;
upd_cm  = 0.*cm;  upd_cx = 0.*cx;
cm = cmq;  cx = cxq;
while rnorm>tol && it<mxit

    Kx = cx./(cm+eps);

    cmK = c    ./(m + x.*Kx + eps); 
    cxK = c.*Kx./(m + x.*Kx + eps); 

    res_cm = cm - cmK./sum(cmK,3);
    res_cx = cx - cxK./sum(cxK,3);

    upd_cm = - 0.95.*res_cm + 0.25.*upd_cm;
    upd_cx = - 0.95.*res_cx + 0.25.*upd_cx;

    cm = max(0,cm + upd_cm);
    cx = max(0,cx + upd_cx);

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

end

% record timing
TCtime = TCtime + toc - eqtime.*Rcouple;
