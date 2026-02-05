% get residual of thermochemical equations
resnorm    = norm(res_S  (:)*dt)./(norm(S  (:))+eps) ...
           + norm(res_C  (:)*dt)./(norm(rho(:))+eps) ...
           + norm(res_M  (:)*dt)./(norm(rho(:))+eps) ...
           + norm(res_X  (:)*dt)./(norm(rho(:))+eps);% ...
           % + norm(res_MFD(:)*dt)./(norm(rho(:))+eps) ...
           % + norm(res_CMP(:)*dt)./(norm(rho(:))+eps);

if iter==1 || resnorm>resnorm0; resnorm0 = resnorm + 1e-32; end  % reset reference residual

% check for solver divergence or failing
if isnan(resnorm); error('!!! Solver failed with NaN: end run !!!'); end

% report iterations
if     iter >=  0  && iter <  10
    fprintf(1,'    ---  iter =    %d;   abs res = %4.4e;   rel res = %4.4e \n',iter,resnorm,resnorm/resnorm0);
elseif iter >= 10  && iter < 100
    fprintf(1,'    ---  iter =   %d;   abs res = %4.4e;   rel res = %4.4e \n',iter,resnorm,resnorm/resnorm0);
elseif iter >= 100 && iter < 1000
    fprintf(1,'    ---  iter =  %d;   abs res = %4.4e;   rel res = %4.4e \n',iter,resnorm,resnorm/resnorm0);
end 

% plot convergence of outer iterations
if plot_cv
    figure(100); if iter==1; else; hold on; end
    plot(iterglobal,log10(resnorm),'k.','MarkerSize',15,'LineWidth',1.5); box on; axis tight;
    drawnow;
end
