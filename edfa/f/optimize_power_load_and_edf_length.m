function [E, Signal, exitflag, num, approx] = optimize_power_load_and_edf_length(method, E, Pump, Signal, problem, verbose)
%% Optimize power loading and EDF length
% Inputs
% - method: optimization method. 
% Assuming that channels can be either on or off: method = 'interp' or 'fminbnd'
% Assuming that channels power is upper bounded by Pon: method = 'particle swarm'
% - E: isntance of class EDF
% - Pump: instance of class Channels corresponding to pump
% - Signal: instance of class Channels corresponding to signals
% - problem: struct with problem parameters: Pon, Namp, spanAttdB, and df
% Outputs:
% - Lopt: optimal EDF length
% - Signal: instance of class Channels corresponding to signals.
% - exitflag: reason optimizatoin ended. Meaning depends on each algorithm. 
% exitflag = 1 usually means that relative changed in objective function
% was smaller than tolerance
% - num: struct containing {GaindB, Pase, SE, SNRdB} calculated numerically
% - approx: struct containing {GaindB, Pase, SE, SNRdB} calculated
% analytically 

% Unpack problem parameters
Pon = problem.Pon;
PoffdBm = Watt2dBm(eps); % assign small power to off channels
spanAttdB = problem.spanAttdB;
% Excess noise
problem.excess_noise = E.analytical_excess_noise(Pump, Signal);
if isfield(problem, 'excess_noise_correction')
    problem.excess_noise = problem.excess_noise*problem.excess_noise_correction;
end


%% Optimization
fprintf('Optimization method = %s\n', method)
switch lower(method)
    case 'fminbnd'
        %% Find EDF length that maximizes the number of ON channels
        % This uses Matlab's fminbnd function to search for the optimal fiber length
        % (Doesn't produce very consistent results)
        
        [Lopt, ~, exitflag] = fminbnd(@(L) -sum(max_channels_on(L, E, Pump, Signal, Pon, spanAttdB)), 0, E.maxL);

        if exitflag ~= 1
            warning('optimize_edf_length: optmization exited with exitflag %d\n', exitflag)
        end

        E.L = Lopt;

        % Evaluate objective at optimal EDF length
        onChs = max_channels_on(Lopt, E, Pump, Signal, Pon, spanAttdB);
        Signal.P = zeros(1, Signal.N);
        Signal.P(onChs) = Pon;

    case 'interp'
        %% Find EDF length that maximizes the number of ON channels
        % Optimal EDF length is found by interpolation. Number of ON
        % channels is calculated for "Nsteps" points from EDF.L = 1 to
        % EDF.Lmax. Optimal result is obtained by interpolation
        Nsteps = 40;
        Ledf = linspace(1, E.maxL, Nsteps); % EDF is at least a meter long
        BW = zeros(size(Ledf));
        parfor k = 1:Nsteps
            BW(k) = sum(max_channels_on(Ledf(k), E, Pump, Signal, Pon, spanAttdB));
        end

        Lfit = linspace(1, E.maxL, 1e3);
        [~, idx] = max(spline(Ledf, BW, Lfit));
        Lopt = Lfit(idx);

        if exist('verbose', 'var') && verbose
            figure(202), clf, hold on, box on
            stem(Ledf, BW)
            plot(Lopt, max(BW), 'or')
            xlabel('EDF length (m)')
            ylabel('Numer of channels On')
        end
        
        E.L = Lopt;
        
        % Evaluate objective at optimal EDF length
        onChs = max_channels_on(Lopt, E, Pump, Signal, Pon, spanAttdB);
        Signal.P = zeros(1, Signal.N);
        Signal.P(onChs) = Pon;
        exitflag = 1;
    case 'particle swarm'
        %% Optimize EDF length and power allocation jointly
        if isfield(problem, 'SwarmSize')
            SwarmSize = problem.SwarmSize;
        else
            SwarmSize = min(200, 20*(Signal.N+1));
        end
                            
        la = [0 -Inf*ones(1, Signal.N)]; % lower bound
        lb = [E.maxL Watt2dBm(Pon)*ones(1, Signal.N)]; % upper bound
        options = optimoptions('particleswarm', 'Display', 'iter', 'UseParallel', true,...
                'MaxStallTime', 60, 'MaxStallIterations', 100, 'SwarmSize', SwarmSize,...
                'InitialSwarmSpan', [20, 10*ones(1, Signal.N)]); % SwarmSpan for fiber length is 20m and 30 dB for signal power
                
        if isfield(problem, 'nonlinearity') && problem.nonlinearity
            disp('IMPORTANT: optimization includes fiber nonlinearity')

            [X, relaxed_SE, exitflag] = particleswarm(@(X) capacity_nonlinear_regime_relaxed(X, E, Pump, Signal, problem),...
                Signal.N+1, la, lb, options);
        else
            [X, relaxed_SE, exitflag] = particleswarm(@(X) -capacity_linear_regime_relaxed(X, E, Pump, Signal, problem),...
            Signal.N+1, la, lb, options);
        end
        
        fprintf('- Relaxed objective: Total Capacity = %.3f (bits/s/Hz)\n', -relaxed_SE)
        
        E.L = X(1);
        Signal.P = dBm2Watt(X(2:end));
        GaindB = E.semi_analytical_gain(Pump, Signal);
        Signal.P(GaindB <= spanAttdB) = 0; % turn off channels that don't meet gain requirement
      
    case 'local'
        %% Local optimizaiton using a gradient-based algorithm. Results from particle swarm optmization should be provided as starting point
        options = optimoptions('fmincon', 'Algorithm', 'trust-region-reflective',...
            'Display', 'iter', 'UseParallel', true,...
            'CheckGradients', true, 'SpecifyObjectiveGradient', true, 'FiniteDifferenceType', 'central',...
            'MaxFunctionEvaluations', 1e4);
        
        Signal.P(Signal.P == 0) = eps; % to respect the bounds
        la = [0 PoffdBm*ones(1, Signal.N)]; % lower bound
        lb = [E.maxL Watt2dBm(Pon)*ones(1, Signal.N)]; % upper bound
        [X, ~, exitflag] = fmincon(@(X) capacity_nonlinear_regime_relaxed(X, E, Pump, Signal, problem), ...
            [E.L Signal.PdBm], [], [], [], [], la, lb, [], options);

        E.L = X(1);
        Signal.P = dBm2Watt(X(2:end));
        GaindB = E.semi_analytical_gain(Pump, Signal);
        Signal.P(GaindB <= spanAttdB) = 0; % turn off channels that don't meet gain requirement
    case 'none' % doesn't do anything. Just use this function for plotting
        exitflag = 1;
        GaindB = E.semi_analytical_gain(Pump, Signal);
        Signal.P(GaindB <= spanAttdB) = 0; % turn off channels that don't meet gain requirement
    otherwise
        error('optimize_power_load_and_edf_length: invalid method')
end

fprintf('- Optimal EDF length = %.2f\n', E.L)
NChOn = sum(Signal.P ~= 0);
fprintf('- Number of channels ON = %d\n', NChOn)

% 
offChs = (Signal.P == 0);
Signal.P(offChs) = eps; % set to small power to calculate gain

% Compute capacity using numerical and semi-analytical (approx) methods 
if isfield(problem, 'nonlinearity') && problem.nonlinearity
    [num, approx] = capacity_nonlinear_regime(E, Pump, Signal, problem);
    [~, ~, SElamb_relaxed] = capacity_nonlinear_regime_relaxed([E.L Signal.P], E, Pump, Signal, problem);
else
    [num, approx] = capacity_linear_regime(E, Pump, Signal, problem);
    [~, SElamb_relaxed] = capacity_linear_regime_relaxed([E.L Signal.P], E, Pump, Signal, problem);
end

Signal.P(offChs) = 0;

% 
fprintf('- Numerical: Total Capacity = %.3f (bits/s/Hz) | Avg. Capacity = %.3f (bits/s/Hz)\n',...
    sum(num.SE), sum(num.SE)/NChOn)
fprintf('- Approximated: Total Capacity = %.3f (bits/s/Hz) | Avg. Capacity = %.3f (bits/s/Hz)\n',...
    sum(approx.SE), sum(approx.SE)/NChOn)

SignalOut = Signal;
SignalOut.P = dBm2Watt(Signal.PdBm + num.GaindB);
[PCE, PCEmax] = E.power_conversion_efficiency(Pump, Signal, SignalOut);
fprintf('- Power conversion efficiency = %.2f %% (max = %.2f %%)\n', 100*PCE, 100*PCEmax)
        
% Plot    
if exist('verbose', 'var') && verbose   
    if length(spanAttdB) == 1
        spanAttdB = spanAttdB*ones(size(Signal.wavelength));
    end
    
    % Plot Results
    lnm = Signal.wavelength*1e9;
    figure(203), hold on, box on
    hplot = plot(lnm, approx.GaindB, 'DisplayName', 'Approximated');
    plot(lnm, num.GaindB, '--', 'Color', get(hplot, 'Color'),'DisplayName', 'Numerical')
    plot(lnm(not(offChs)), approx.GaindB(not(offChs)), 'ob', 'DisplayName', 'Channels on')
    plot(lnm(offChs), approx.GaindB(offChs), 'xr', 'DisplayName', 'Channels off')
    plot(lnm, spanAttdB, 'k', 'DisplayName', 'Span attenuation')
    axis([lnm([1 end]), 0, 20])
    xlabel('Wavelength (nm)')
    ylabel('Gain (dB)')
    if isempty(legend(gca))
    	legend('-DynamicLegend', 'Location', 'SouthEast')
    end
    title(sprintf('EDF %s, L = %.2f m', E.type, E.L), 'Interpreter', 'none')
    
    figure(204), hold on, box on
    hplot = plot(lnm, Watt2dBm(approx.Pase), 'DisplayName', 'Approximated');
    plot(lnm, Watt2dBm(num.Pase), '--', 'Color', get(hplot, 'Color'), 'DisplayName', 'Numerical')
    xlabel('Wavelength (nm)')
    ylabel('ASE (dBm)')
    if isempty(legend(gca))
        legend('-DynamicLegend', 'Location', 'SouthEast')
    end
    title(sprintf('EDF %s, L = %.2f m', E.type, E.L), 'Interpreter', 'none')
    
    figure(205), hold on, box on
    hplot = plot(lnm, approx.SNRdB, 'DisplayName', 'Approx');
    plot(lnm, num.SNRdB, '--', 'Color', get(hplot, 'Color'), 'DisplayName', 'Numerical')
    xlabel('Wavelength (nm)')
    ylabel('SNR (dB)')
    legend('-DynamicLegend', 'Location', 'SouthEast')
    title(sprintf('EDF %s, L = %.2f m', E.type, E.L), 'Interpreter', 'none')
    axis([lnm(1) lnm(end) 0 25])
    
    figure(206), hold on, box on
    hplot = plot(lnm, approx.SE, 'DisplayName', 'Approximated');
    plot(lnm, SElamb_relaxed, ':', 'Color', get(hplot, 'Color'), 'DisplayName', 'Relaxed')
    plot(lnm, num.SE, '--', 'Color', get(hplot, 'Color'), 'DisplayName', 'Numerical')
    xlabel('Wavelength (nm)')
    ylabel('Capacity (bits/s/Hz)')
    if isempty(legend(gca))
        legend('-DynamicLegend', 'Location', 'SouthEast')
    end
    title(sprintf('EDF %s, L = %.2f m', E.type, E.L), 'Interpreter', 'none')
       
    figure(207), 
    subplot(211), box on, hold on
    plot(lnm, Signal.PdBm)
    xlabel('Wavelength (nm)')
    ylabel('Power loading (dBm)')

    subplot(212), box on, hold on
    plot(lnm, Signal.PdBm + approx.GaindB)
    xlabel('Wavelength (nm)')
    ylabel('Output power (dBm)')
    
    if isfield(problem, 'nonlinearity') && problem.nonlinearity
        figure(208), hold on, box on
        plot(lnm, Watt2dBm(num.NL), 'DisplayName', 'Nonlinear noise')
        plot(lnm, Watt2dBm(num.Pase), 'DisplayName', 'ASE')
        xlabel('Wavelength (nm)')
        ylabel('Noise power (dBm)')
        legend('-DynamicLegend', 'Location', 'SouthEast')
    end
    drawnow
end