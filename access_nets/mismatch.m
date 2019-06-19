%% 25 km CD precompensation with 10 km fiber

load results/preeq.mat optgain

if exist('results/mismatch.mat','file')
    load results/mismatch
else
    %% simulate with optimum APD gain
    PtxdBm = -35:-25;
    Poff = [0 5 6];
    L = [10 0];

    for k = length(L):-1:1
        for j=length(M):-1:1
            ber(k,j) = pam_access_nets_qsub(30,M(j),L(k),1358,30,'MZM','optimized',2,optgain(j,1),22,'precomp25',1e-3,PtxdBm+Poff(j),0.8,3);
        end
    end
    %% receiver sensitivities

    BERtarget = [0.00145 0.0032 0.0052 0.0067 0.008 0.0095 0.012 0.015 0.016 0.019 0.032];
    BERtarget = BERtarget(end:-1:1);

    Preqenum = zeros(length(BERtarget),length(M),length(L));
    Preq = zeros(length(BERtarget),length(M),length(L));

    for k = 1:length(L)
        for im=1:length(M)
            m = M(im);
                BERcount = ber(k,im).count;
                BERenum = ber(k,im).enum;
                PrxdBm = ber(k,im).PrxdBm;

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

    % save('results/mismatch.mat','M','L','ber','BERtarget','Preqenum','Preq','PtxdBm','Poff')
end
%%
figure(3)
preeq = load('results/preeq.mat','Preq');
h0=plot(log10(BERtarget),preeq.Preq(:,:,1),'+-');
hold on
set(gca,'colororderindex',1)
h1=plot(log10(BERtarget),preeq.Preq(:,:,2),'.');
set(gca,'colororderindex',1)
h2=plot(log10(BERtarget),Preq(:,:,1),'+:');
legend([h0; h1(1); h2(1)],'OOK 25 km','2DPS 3-PAM','4-PAM', 'b2b','10km mismatch')
grid on
ylabel('P_{rec,req} (dBm)')
xlabel('log_{10} BER')
hold off
