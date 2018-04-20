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

%% ============================ Receiver ==================================
% APD
% apd(GaindB, ka, BW, R, Id)
APD = apd(11, 0.5, Inf, 1, 10e-9);

%% TIA-AGC
% One-sided thermal noise PSD
N0 = (30e-12).^2; 

%% find optimal APD gain, with corresponding thermal and shot noise

G = linspace(5,50);
thermal = sqrt(N0*Rs);
aopt = cell(3,1);
bopt = cell(3,1);
PrecreqdB = zeros(3);

j=1;
figure(1)
for M=[2 4 8]
    levels = zeros(length(G),M);
    thresh = zeros(length(G),M-1);
    shot = zeros(length(G),M);
    Prec = zeros(length(G),1);
%     Prec_2pam = zeros(length(G),1);

    Pe = log2(M)*BERtarget*M/(2*(M-1));
    Qreq = fzero(@(Q) qfunc(Q) - Pe, 0);
    for i=1:length(G)
        APD.Gain = G(i);
        noise_std =@(P) sqrt(thermal^2+APD.varShot(P/APD.Geff,Rs));

        for level = 1:M-1
            % Find threshold
            thresh(i,level) = levels(i,level) + abs(Qreq*noise_std(levels(i,level)));

            % Find next level  
            [dPlevel, ~, exitflag] = fzero(@(dPlevel) qfunc(abs(dPlevel)/noise_std(thresh(i,level) + abs(dPlevel))) - Pe, 0);    

            if exitflag ~= 1
                warning('level_spacing_optm: level optimization did not converge');     
            end

            levels(i,level+1) = thresh(i,level) + abs(dPlevel);
        end
        Prec(i) = 10*log10(mean(levels(i,:))/APD.Geff*1e3);
    %     Prec_2pam(i) = 10*log10(Qreq/APD.R*(APD.q*APD.Fa*Qreq*Rs+thermal/APD.Gain)*1e3);
        shot(i,:) = sqrt(APD.varShot(levels(i,:)/APD.Geff,Rs));
    end
    [PrecreqdB(j),imin] = min(Prec);
    aopt{j} = levels(imin,:);
    bopt{j} = thresh(imin,:);
    
    yyaxis right
    h = plot(G,shot(:,end)/thermal);
    hasbehavior(h,'legend',false)
    yyaxis left
    plot(G,Prec,'displayname',sprintf('%dPAM',M))
    hold on
    h = plot(G(imin),PrecreqdB(j),'x');
    hasbehavior(h,'legend',false)
%     plot(ones(2,1)*G(Prec==min(Prec)),ylim)
%     hold off
    xlabel('APD gain')
    j = j+1;
end
hold off
yyaxis left
ylabel('Receiver sensitivity (dBm)')
yyaxis right
ylabel('shot/thermal noise ratio')
ylim([0 10])
legend('-dynamiclegend')
%% 2-PAM equations
% ka = 0.1:0.02:1;
% Mopt = sqrt((thermal./(Qreq*APD.q*Rs)+ka-1)./ka);
% FA = ka.*Mopt+(1-ka).*(2-1./Mopt);
% Precreq = 2*APD.q*Rs/APD.R*Qreq^2*(ka.*Mopt+1-ka);
% PrecdB2 = 10*log10(Qreq/APD.R*(APD.q*FA*Qreq*Rs+thermal./Mopt)*1e3);
% PrecdB = 10*log10(Precreq*1e3);
% shot_std = sqrt(2*APD.q*Mopt.^2.*FA*2*APD.R.*Precreq*Rs);
% 
% figure(2)
% yyaxis left
% plot(ka, Mopt)
% ylabel('optimal gain')
% yyaxis right
% plot(ka, shot_std/thermal)
% ylabel('shot/thermal noise ratio')
% % plot(ka, [PrecdB; PrecdB2]')
% % ylabel('Receiver sensitivity (dBm)')
% hold on
% ybounds = ylim;
% plot(ones(2,1)*APD.ka,ybounds,'k')
% ylim(ybounds)
% hold off
% xlabel('k_A')
% title('2PAM theoretical')
