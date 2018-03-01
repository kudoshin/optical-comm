function [G, Prec, thermal, shot, levels, thresh] = apd_gain(M)
%% Calculate BER of IM-DD system using M-PAM
% - Equalization is done using a fractionally spaced linear equalizer
% Simulations include modulator, fiber, optical amplifier (optional) characterized 
% only by gain and noise figure, optical bandpass filter, antialiasing 
% filter, sampling, and linear equalization

addpath ../f  % general functions
addpath ../apd % for PIN photodetectors
addpath ../mpam

%% Simulation parameters
BERtarget = 1.8e-4; 
Rs = 12.5e9;

%% Pulse shape
Tx.pulse_shape = select_pulse_shape('rect', 1);

%% M-PAM
% PAM(M, bit rate, leve spacing : {'equally-spaced', 'optimized'}, pulse
% shape: struct containing properties of pulse shape 
mpam = PAM(M, Rs*log2(M), 'equally-spaced', Tx.pulse_shape);

%% ============================ Receiver ==================================
% APD
% apd(GaindB, ka, BW, R, Id)
APD = apd(11, 0.5, Inf, 1, 10e-9);

%% TIA-AGC
% One-sided thermal noise PSD
N0 = (30e-12).^2; 

%% find optimal APD gain, with corresponding thermal and shot noise
%
Pe = log2(M)*BERtarget*M/(2*(M-1));
G = linspace(1,40,10);
thermal = sqrt(N0*Rs);
levels = zeros(length(G),M);
thresh = zeros(length(G),M-1);
shot = zeros(length(G),M);
Prec = zeros(length(G),1);
for i=1:length(G)
    APD.Gain = G(i);
    noise_std =@(P) sqrt(thermal^2+APD.varShot(P/APD.Geff,Rs));
    
    mpam = mpam.optimize_level_spacing_gauss_approx(BERtarget, -Inf, noise_std, true);
    for level = 1:M-1
        % Find threshold
        sig = noise_std(levels(i,level));

        [dPthresh, ~, exitflag] = fzero(@(dPthresh) qfunc(abs(dPthresh)/sig) - Pe, 0);

        if exitflag ~= 1
            warning('level_spacing_optm: threshold optimization did not converge');
        end

        thresh(i,level) = levels(i,level) + abs(dPthresh);

        % Find next level  
        [dPlevel, ~, exitflag] = fzero(@(dPlevel) qfunc(abs(dPlevel)/noise_std(thresh(i,level) + abs(dPlevel))) - Pe, 0);    

        if exitflag ~= 1
            warning('level_spacing_optm: level optimization did not converge');     
        end

        levels(i,level+1) = thresh(i,level) + abs(dPlevel);
    end
    Prec(i) = 10*log10(mean(levels(i,:))/APD.Geff*1e3);
    shot(i,:) = sqrt(APD.varShot(levels(i,:)/APD.Geff,Rs));
end

%%
figure(1)
yyaxis left
plot(G,Prec)
ylabel('Receiver sensitivity (dBm)')
yyaxis right
semilogy(G,shot/thermal)
ylabel('noise std dev')
xlabel('APD gain')
