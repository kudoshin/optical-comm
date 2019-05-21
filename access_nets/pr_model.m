PrecdB = -30:1.25:-10;
Prec = 10.^(PrecdB/10)*1e-3;
pdf = [0.02 0.025 0.04 0.05 0.12 0.2 0.475 0.72 1 0.89 0.65 0.34 0.2 0.09 0.025 0 0];
s = 0.9e-5;
rayleigh = Prec.*exp(-Prec.^2/(2*s^2))/s^2;
Pav = -20;
sig = 2.5;
gaussian = 1/sqrt(2*pi*sig^2)*exp(-(PrecdB-Pav).^2/(2*sig^2));
pdf = pdf*max(gaussian);
% plot(Prec,[pdf;rayleigh;gaussian]')
% L = 0:50;
% s = 20;
% rayleigh = L.*exp(-L.^2/(2*s^2))/s^2;
% plot(L,rayleigh)
% 
% loss = 0.3;
% Prec = -10-loss*L;
% plot(Prec,rayleigh)
plot(PrecdB,[pdf;rayleigh.*Prec*log(10)/10;gaussian]')

%% splice loss + distance model
dt = 0.0002;
tdB = -10:dt:0;
t = 10.^(tdB/10);
alpha = 0.096;
beta = 0.265;
psplice = zeros(size(t));
r1 = -log(t)<alpha^2;
r2 = -log(t)<beta^2;
r3 = -log(t)<alpha^2+beta^2;
r4 = r2 & ~r1;
r5 = r3 & ~r2;
psplice(r1) = pi./(4*t(r1)*alpha*beta);
psplice(r4) = 1./(2*t(r4)*alpha*beta).*atan(alpha./sqrt(-log(t(r4))-alpha^2));
psplice(r5) = 1./(2*t(r5)*alpha*beta).*(atan(alpha./sqrt(-log(t(r5))-alpha^2))-acot(beta./sqrt(-log(t(r5))-beta^2)));
% plot(t,psplice)
pdBsplice = psplice.*t*log(10)/10;
Nsplice = 3;
% plot(10*log10(t),pdBsplice)
Nconn1 = 3;
Nconn2 = 1;
meangauss = -Nconn1*0.3+Nconn2*0.5;
pdBgauss = normpdf(-3:dt:2,meangauss,0.5);

loss = 0.4; % dB/km
meandist = 8; % km
meanloss = meandist*loss; % dB
s2 = meanloss^2*2/pi; %km^2
pdistdB = -tdB/s2.*exp(-tdB.^2/(2*s2));
% plot(tdB,pdistdB)

pdfVar = pdBgauss;
for i=1:Nsplice
    pdfVar = conv(pdBsplice,pdfVar)*dt;
end
pdfTot = conv(pdfVar,pdistdB)*dt;
cdfTot = cumsum(pdfTot)/sum(pdfTot);
%%
lossdB = -12:dt:2;
pdfVar = pdfVar(end-length(lossdB)+1:end);
pdfTot = pdfTot(end-length(lossdB)+1:end);
cdfTot = cdfTot(end-length(lossdB)+1:end);
%%
fig = figure(1);
fig.Position = [0 550 560 260];
yyaxis left
% set(0,'defaultaxesfontsize',16)
% set(0,'defaultaxesfontname','times new roman')
h=plot(lossdB,pdfTot,'linewidth',1.5,'color','k');
% hold on
yyaxis right
plot(lossdB,cdfTot,'linewidth',2)
hold off
xlabel('Relative loss (dBm)')
xlim([-12 0])
% ylim([0 1])
grid on
legend('pdf','cdf')
ylim([0 1])
%%
save prmodel lossdB pdfTot cdfTot pdfVar