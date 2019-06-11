%% 25 km, b2b
if exist('results\posteq.mat','file')
    load results\posteq.mat
else
    apdgain = 12.2:0.1:17;
    M = [2 3 4];
    L = [25 0];
    [PreqdBm,ber20dbm] = deal(zeros(length(apdgain),length(M),length(L)));
    for k=1:length(L)
        for i=1:length(apdgain)
            for j=1:length(M)
                ber = pam_access_nets_qsub(30,M(j),L(k),1358,30,'MZM','optimized',2,apdgain(i),22,'fixed td-sr-le',1e-3,-20,0.8,0);
                PreqdBm(i,j,k) = 10*log10(ber.preq*1e3);
                ber20dbm(i,j,k) = ber.enum(end);
            end
        end
    end
    [minP, imin] = min(PreqdBm,[],1);
    minP = squeeze(minP);
    optgain = squeeze(apdgain(imin));
    save('results/posteq.mat','M','L','apdgain','PreqdBm','ber20dbm','optgain','minP')
end
% 
%%
figure(1)
plot(apdgain,PreqdBm(:,:,1))
hold on
set(gca,'colororderindex',1)    
plot(apdgain,PreqdBm(:,:,2),'--')
hold off
grid on 
xlabel('APD gain (dB)')
ylabel('Receiver sensitivity (dBm)')
legend('2-PAM 25km', '2DPS 3-PAM', '4-PAM','b2b','location','northwest')

%% simulate with optimum APD gain
PtxdBm = -35:-25;
Poff = [0 4 6];
for k=length(L):-1:1
    for j=length(M):-1:1
        ber(j,k) = pam_access_nets_qsub(30,M(j),L(k),1358,30,'MZM','optimized',2,optgain(j,k),22,'fixed td-sr-le',1e-3,PtxdBm+Poff(j),0.8,2);
    end
end

%% receiver sensitivities

BERtarget = [0.00145 0.0032 0.0052 0.0067 0.008 0.0095 0.012 0.015 0.016 0.019 0.032];
BERtarget = BERtarget(end:-1:1);
        
Preqenum = zeros(length(BERtarget),length(M),length(L));
Preq = zeros(length(BERtarget),length(M),length(L));

for im=1:length(M)
    m = M(im);
    for k=1:length(L) 
            BERcount = ber(im,k).count;
            BERenum = ber(im,k).enum;
            PrxdBm = ber(im,k).PrxdBm;
            
            idx = find(BERenum <= 0.1 & BERenum >= 1e-4);
            for j=1:length(BERtarget)
                [Preqenum(j,im,k), f] = fit_ber(PrxdBm(idx), BERenum(idx), BERtarget(j));
                Preq(j,im,k) = fit_ber(PrxdBm(idx), BERenum(idx), BERtarget(j));
            end
            % Plot
            figure(5), clf, hold on, box on
            h = plot(PrxdBm, f(PrxdBm), '-');
            plot(PrxdBm, log10(BERcount), '--', 'Color', get(h, 'Color'));
            plot(PrxdBm, log10(BERenum), '-.', 'Color', get(h, 'Color'));
            xlim(PrxdBm([1 end]))
            ylim([-5 0])
            grid on
            hline(log10(BERtarget),'r--')
            hold on
            drawnow  
    end
end
hold off

% save('results/posteq.mat','M','L','apdgain','PreqdBm','ber20dbm','optgain','minP','ber','BERtarget','Preqenum','Preq','PtxdBm','Poff')
%%
figure(3)
plot(log10(BERtarget),Preq(:,:,1),'x-')
hold on
set(gca,'colororderindex',1)
plot(log10(BERtarget),Preq(:,:,2),'x--')
legend('OOK 25 km','2DPS 3-PAM','4-PAM', 'b2b')
grid on
ylabel('P_{rec,req} (dBm)')
xlabel('log_{10} BER')
hold off
