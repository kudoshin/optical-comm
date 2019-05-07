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
tdB = -6:dt:0;
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
% plot(10*log10(t),pdBsplice)
pdBgauss = normpdf(-2:dt:2,0,0.5);

meandist = 4; % km
s2 = meandist^2*2/pi; %km^2
pdistdB = -4*tdB/s2.*exp(-4*tdB.^2/(2*s2));
% plot(tdB,pdistdB)

pdfTot = conv(conv(conv(pdBsplice,pdBsplice),pdBsplice),conv(pdBgauss,pdistdB))*dt^4;
cdfTot = cumsum(pdfTot)/sum(pdfTot);
lossdB = -12:dt:0;
pdfTot = pdfTot(end-length(lossdB)+1:end);
cdfTot = cdfTot(end-length(lossdB)+1:end);
%%
fig = figure(1);
fig.Position = [0 550 560 260];
set(0,'defaultaxesfontsize',16)
set(0,'defaultaxesfontname','times new roman')
yyaxis left
h=plot(lossdB,pdfTot,'linewidth',1.5,'color','k');
ylabel('pdf')
yyaxis right
plot(lossdB,cumsum(pdfTot)*dt,'linewidth',2)
ylabel('cdf')
xlabel('Relative loss (dBm)')
xlim([-10 0])
% ylim([0 1])
grid on
% save prmodel lossdB pdfTot cdfTot