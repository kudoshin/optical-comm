posteq = load('posteq.mat');
preeq = load('preeq.mat');
mismatch = load('mismatch.mat');
noeq = load('noeq.mat');

%%
figure(1)
h0=semilogy(1,1,1,1,1,1,1,1,'k',1,1,'k--',1,1,'k:',0,0,'k+',0,0,'kx');
hold on
set(gca,'colororderindex',1)
semilogy(posteq.PtxdBm'+posteq.Poff,[posteq.ber(1,:).enum],'x-')
set(gca,'colororderindex',1)
semilogy(posteq.PtxdBm'+posteq.Poff,[posteq.ber(2,:).enum],'x--')
% set(gca,'colororderindex',1)
% semilogy(posteq.PtxdBm'+posteq.Poff,[posteq.ber(3,:).enum],'x:')
set(gca,'colororderindex',1)
semilogy(preeq.PtxdBm'+preeq.Poff,[preeq.ber(1,:).enum],'+-')
% set(gca,'colororderindex',1)
% semilogy(preeq.PtxdBm'+preeq.Poff,[preeq.ber(2,:).enum],'+--')
set(gca,'colororderindex',1)
semilogy(mismatch.PtxdBm'+mismatch.Poff,[mismatch.ber(1,:).enum],'+:')
set(gca,'colororderindex',1)
semilogy(noeq.PtxdBm'+noeq.Poff,[noeq.ber(1,:).enum])
set(gca,'colororderindex',1)
semilogy(noeq.PtxdBm'+noeq.Poff,[noeq.ber(2,:).enum],'--')
hold off
grid on
ylim(10.^[-4 -1])
xlim([-32 -19])
hline(noeq.BERtarget([2 end]))
xlabel('Received optical power (dBm)')
ylabel('Bit error rate')

legend(h0,'OOK','3-PAM','4-PAM','25 km','b2b','10 km','Tx eq.','Rx eq.','location','best')

