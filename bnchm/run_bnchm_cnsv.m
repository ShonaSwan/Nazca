% prepare workspace
clear; close all;

% load default parameters
run('../usr/par_default')

% test decreasing time step
ATOL = [1e-4,1e-7,1e-10];

for atol = ATOL

    % set run parameters
    runID    =  'bnchm_cnsv';        % run identifier
    nop      =  20;                  % output frame plotted/saved every 'nop' time steps
    plot_op  =  1;                   % switch on to live plot of results
    plot_cv  =  0;                   % switch on to live plot iterative convergence
    save_op  =  0;

     % set model domain parameters
    D        =  200e3;                  % chamber depth [m]
    L        =  1*D;                  % chamber width [m]
    N        =  100;                 % number of grid points in z-direction (incl. 2 ghosts)
    h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
    alpha    =  0.4;                % iterative step size parameter

    % set model timing parameters
    Nt       =  nop;                 % number of time steps to take
    dt       =  10*yr;                   % set initial time step
    maxit    =  100;                  % maximum outer its
    minit    =  0.05;                % maximum initial melt fraction (Initial reduction of melt)
    mumax    =  0.20;                 % Setting upper limit for melt fraction in 

    % model set up switches (plume or MOR)
    init_mode   =  'plume';               % 'plume' or 'MOR'
    bndmode     =  1;                   % boundary assimilation mode (0 = MOR; 1 = Plume 
    meansw      =  0;                   % 0 = Geometric mean 1 = Arithmetic mean
    erupt_ratio = 0.5;                % 1 = all eruption (surface), 0 = all emplacement (intrusion at moho), values in between = partitioning

    % set initial thermo-chemical state of the Plume 
    dT_plume  = 200;                                % Temperature difference between the plume and the mantle 
    pl_width  = 50e3;                               % Width of the plume [m]
    pl_local  = L/2; % L/2 + 100                    % Location of the mantle plume along the bottom boundary [m]
    c_plume   = [0.80 0.18 0.02 0];                 % components of plume (maj comp, H2O) [wt] (will be normalised to unit sum!)
    trc_plume = [10.0, 10.0, 2.0, 0.1, 0.1, 2.0];   % trace elements system plume [wt ppm]

    % create output directory
    if ~isfolder([outdir,'/',runID])
        mkdir([outdir,'/',runID]);
    end

    % run code
    run('../src/main')

    % plot convergence
    ES = norm(diff(hist.ES(Nt/2:Nt  ))./diff(hist.time(Nt/2:Nt)),'fro')./sqrt(length(diff(hist.ES(Nt/2:Nt))         ));
    EB = norm(diff(hist.EB(Nt/2:Nt  ))./diff(hist.time(Nt/2:Nt)),'fro')./sqrt(length(diff(hist.EB(Nt/2:Nt))         ));
    EM = norm(diff(hist.EM(Nt/2:Nt  ))./diff(hist.time(Nt/2:Nt)),'fro')./sqrt(length(diff(hist.EM(Nt/2:Nt))         ));
    EX = norm(diff(hist.EX(Nt/2:Nt  ))./diff(hist.time(Nt/2:Nt)),'fro')./sqrt(length(diff(hist.EX(Nt/2:Nt))         ));
    EC = norm(diff(hist.EC(Nt/2:Nt,:))./repmat(diff(hist.time(Nt/2:Nt)),cal.ncmp,1).','fro')./sqrt(length(diff(hist.EC(Nt/2:Nt))*cal.ncmp));
    ET = norm(diff(hist.ET(Nt/2:Nt,:))./repmat(diff(hist.time(Nt/2:Nt)),cal.ntrc,1).','fro')./sqrt(length(diff(hist.ET(Nt/2:Nt))*cal.ntrc));

    clist = [colororder;[0 0 0]];

    fh20 = figure(20);
    subplot(4,1,1);
    semilogy((hist.time(Nt/2:Nt-1)+hist.time(Nt/2+1:Nt))/2,abs(diff(hist.ES(Nt/2:Nt)))./diff(hist.time(Nt/2:Nt)),'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,2);
    semilogy((hist.time(Nt/2:Nt-1)+hist.time(Nt/2+1:Nt))/2,abs(diff(hist.EB(Nt/2:Nt)))./diff(hist.time(Nt/2:Nt)),'-',LW{:}); hold on; axis tight; box on; hold on
    semilogy((hist.time(Nt/2:Nt-1)+hist.time(Nt/2+1:Nt))/2,abs(diff(hist.EM(Nt/2:Nt)))./diff(hist.time(Nt/2:Nt)),'-',LW{:}); hold on; axis tight; box on;
    semilogy((hist.time(Nt/2:Nt-1)+hist.time(Nt/2+1:Nt))/2,abs(diff(hist.EX(Nt/2:Nt)))./diff(hist.time(Nt/2:Nt)),'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_{\bar{\rho}}, \Delta E_{F^i}$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,3);
    semilogy((hist.time(Nt/2:Nt-1)+hist.time(Nt/2+1:Nt))/2,abs(diff(hist.EC(Nt/2:Nt,:)))./repmat(diff(hist.time(Nt/2:Nt)),cal.ncmp,1).','-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_{C_j}$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,4);
    semilogy((hist.time(Nt/2:Nt-1)+hist.time(Nt/2+1:Nt))/2,abs(diff(hist.ET(Nt/2:Nt,:)))./repmat(diff(hist.time(Nt/2:Nt)),cal.ntrc,1).','-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_{\Theta_k}$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    xlabel('Time [s]',TX{:},FS{:});

    fh21 = figure(21);
    p1 = loglog(atol,ES,'+','Color',clist(1,:),'MarkerSize',10,'LineWidth',2); hold on; box on;
    p2 = loglog(atol,EB,'s','Color',clist(2,:),'MarkerSize',10,'LineWidth',2);
    p3 = loglog(atol,EM,'o','Color',clist(3,:),'MarkerSize',10,'LineWidth',2);
    p4 = loglog(atol,EX,'d','Color',clist(4,:),'MarkerSize',10,'LineWidth',2);
    p5 = loglog(atol,EC,'^','Color',clist(5,:),'MarkerSize',10,'LineWidth',2);
    p6 = loglog(atol,ET,'v','Color',clist(6,:),'MarkerSize',10,'LineWidth',2);
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('abs. residual tolerance [1]','Interpreter','latex','FontSize',16)
    ylabel('rel. conservation error rate [1/s]','Interpreter','latex','FontSize',16)
    title('Global conservation with nonlinear convergence','Interpreter','latex','FontSize',20)

    if atol == ATOL(end)
        p7 = loglog(ATOL,eps.*ones(size(ATOL)),'k-' ,'LineWidth',2);  % plot trend for comparison
        legend([p1,p2,p3,p4,p5,p6,p7],{'error $S$','error $\bar{\rho}$','error $M$','error $X$','error $C_j$','error $\Theta_k$','machine prec.'},'Interpreter','latex','box','on','location','southeast')
    end
    drawnow;

end

name = [outdir,'/',runID,'/',runID,'_',TINT,'_',ADVN];
print(fh21,name,'-dpng','-r300','-vector');