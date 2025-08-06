% initialise model run
init;

%%%%  Physical Time Stepping Loop %%%%
while time <= tend && step <= Nt
    
    % Time step info
    timing;

    % Store previous solution
    store;
    
    %%%% Reset residuals and iteration count
    resnorm  = 1;
    resnorm0 = resnorm;
    iter     = 1;
    % if frst; alpha = alpha/2; beta = beta/2; end

    %%%%% Non-Linear Iteration Loop %%%%
    while resnorm/resnorm0 >= rtol/(1 + frst*100) && resnorm >= atol/(1 + frst*10) && iter <= maxit*(1 + frst)
        
        %%%% solve thermo-chemical equations
        thermochem;

        %%%% update non-linear parameters and auxiliary variables
        update;

        %%%% solve fluid-mechanics equations
        fluidmech;

        %%%% update geochemical evolution
        geochem;

        %%%% report convergence
        report;

        %%%% Increment iteration count 
        iter = iter+1; 
   
    %%%% The end of the Non-Linear Iteration Loop %%%%
    end
    
    %%%% Upsate phase equlibrium 
    phseql; 

    % record model history
    history;

    % print model diagnostics
    diagnose;

    % plot model results
    if ~mod(step,nop); output; end
    
    % increment time/step
    time = time+dt;
    step = step+1;
    frst = 0;
    % if frst; alpha = alpha*2; beta = beta*2; frst=0; end
 
%%%% The end of the Physical Time Stepping Loop %%%% 
end

% save final state of model
output;

diary off
