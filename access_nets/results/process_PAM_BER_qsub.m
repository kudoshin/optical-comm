%% Process data saved by PAM_BER_qsub.m
% clear, clc, close all

addpath ../../mpam
addpath ../../f/
addpath ../../apd/
addpath ../../apd/f/
addpath ../

M = [2 4 8];
% folder = 'lamb=1550nm_Ntaps=15_ENOB=5_ros=1.00\ModBW=25GHz\RecBW=25GHz\amp=apd\';
folder = 'ModBW=25GHz\RecBW=25GHz\amp=apd\';
% folder = 'BW=InfGHz\amp=pin\';

% PAM_BER_L=0.5km_lamb=1380nm_ModBW=30GHz_amplified=0_Ntaps=9_ENOB=5_ros=1.25

BERtarget = 1.05e-2;%[1.8e-4 4.73e-4 7.6e-4 1.05e-3 1.4e-3];
Amplified = {'pin', 'soa', 'apd'};
ModBWGHz = 10;
ENOB = 5;
Ntaps = 15;
lamb = 1550;
ros = 2;

LineStyle = {'-', '--'};
Marker = {'o', 's', 'v', '^','d'};
Color = {[51, 105, 232]/255, [153,153,155]/255, [255,127,0]/255};
rates = [239/255 73/89 38/54 23/39 14/30];
% Lkm = 0:0.5:10;
Lkm = [0 10 20 30 40 50];


% Fiber = fiber();
D = zeros(1, length(Lkm));  
PrxdBm = cell(length(M),1);
SNR = zeros(length(M),length(BERtarget));

for im=1:length(M)
    m = M(im);
    top = sprintf('12.5Gbd\\%dPAM\\',m);
    PrxdBm{im} = zeros(length(BERtarget),length(Lkm));
    for k=1:length(Lkm)  
        filename = [top folder sprintf('PAM_BER_L=%dkm.mat',...
            Lkm(k))];  

        try 
            S = load(filename, '-mat');
            D(k) = 1e6*S.Fibers(1).D(S.Tx.Laser.wavelength)*S.Fibers(1).L/1e3;

%             P2SNR =@(P) S.Rx.PD.stdNoise(S.
            % Realizations were already averaged
            BERcount = log10(S.ber.count);
            BERgauss = log10(S.ber.awgn);
            BERenum = log10(S.ber.enum);
            
            idx = find(BERcount <= -1 & BERcount >= -5);
            for j=1:length(BERtarget)
                [PrxdBm{im}(j, k), f] = fit_ber(S.Tx.PtxdBm(idx), S.ber.count(idx), BERtarget(j));
            end
            
            % Plot
            figure(3), clf, hold on, box on
            hline = plot(S.Tx.PtxdBm, BERcount, '-o');
            plot(S.Tx.PtxdBm, f(S.Tx.PtxdBm), '-', 'Color', get(hline, 'Color'));
            plot(S.Tx.PtxdBm, BERgauss, '--', 'Color', get(hline, 'Color'));
            plot(S.Tx.PtxdBm, BERenum, '-.', 'Color', get(hline, 'Color'));
            axis([S.Tx.PtxdBm([1 end]) -8 0])
            title(sprintf('L = %.1f km, D = %.2f', S.Fibers(1).L/1e3, D(k)))
            drawnow   
%             pause
        catch e
            filename
            warning(e.message)
            PrxdBm{im}(j, k) = NaN;
        end
    end
    Prec = 10.^(PrxdBm{im}(:,1)/10)*1e-3;
    % no preamp SNR
    SNR(im,:) = 10*log10(S.Rx.PD.Geff^2*Prec.^2./(S.Rx.N0*S.mpam.Rs+S.Rx.PD.varShot(Prec,S.mpam.Rs))); % without noise enhancement
    figure(1), hold on, box on
    for j=1:length(BERtarget)
        plot(Lkm, reshape(PrxdBm{im}(1,:),1,[]), ['-' Marker{j}], 'displayname', sprintf('%dPAM,R=%.2f', m,rates(j)),'Color',Color{im})
    end
    figure(2), hold on
    for j=1:length(BERtarget)
%         plot(SNR(im,j), S.sim.Rb*rates(j)*1e-9, Marker{j}, 'Color',Color{im})
        h=plot(PrxdBm{im}(j,1), S.sim.Rb*rates(j)*1e-9, Marker{j}, 'color',Color{im},'markerfaceColor',Color{im},'markersize',15);
        hasbehavior(h,'legend',false)
    end
%     plot(SNR(im,:),S.sim.Rb*rates*1e-9,'color',Color{im})
    plot(PrxdBm{im}(:,1),S.sim.Rb*rates*1e-9,'color',Color{im},'displayname',sprintf('%d-PAM',m),'linewidth',2)
end
set(gca, 'fontsize', 20)

figure(1)
% xlabel('Dispersion (ps/nm)')
xlabel('Transmission Distance (km)')
ylabel('Receiver sensitivity (dBm)')
legend('-dynamiclegend')
ylim([-22 0])
title('PIN')

figure(2)
% xlabel('SNR')
xlabel('P_{rec,req} (dBm)','fontsize',24)
ylabel('Rate (Gb/s)','fontsize',24)
lgd = legend('-dynamiclegend');
legend('location','northwest')
lgd.FontSize = 20;

figure(2)
% SNR = cellfun(@(c)S.Rx.PD.Geff^2*c(:,1).^2,PrxdBm,'uniformoutput',false)
