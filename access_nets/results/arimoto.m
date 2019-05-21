q = 1.602176e-19;

GdB = 18;
G = 10^(GdB/10);
R = 0.7;
kA = 0.15;
FA = kA*G+(1-kA)*(2-1/G);
N0 = 30e-12^2;

Rs = 25e9;
df = Rs/2;  % half of symbol rate
% df = Rs;

%       _______
% j -> |channel| -> i

PdBm = -50:-10;
C = zeros(1,length(PdBm));

m = 2000;
n = 2000;

power = zeros(n,length(PdBm));
px = zeros(n,length(PdBm));

parfor k = 1:length(PdBm)
    preq = 10^(PdBm(k)/10)*1e-3;
    j = linspace(0,preq*20,n);
    power(:,k) = j;

    nth = N0*df;
    nsh = 2*q*G^2*R*FA*j*df;
    nstd = sqrt(nth+nsh);

    i = linspace(-5*sqrt(nth),G*R*j(end)+5*nstd(end),m)';
    edgesL = [-inf; (i(2:end)+i(1:end-1))/2];
    edgesR = [(i(2:end)+i(1:end-1))/2; inf];

    P = normcdf((edgesR-G*R*j)./nstd) - normcdf((edgesL-G*R*j)./nstd);
    Pu = normcdf((edgesL-G*R*j)./nstd,'upper') - normcdf((edgesR-G*R*j)./nstd,'upper');
    P(P==0)=Pu(P==0);

    H =@(p) -nansum(p.*log(p));

    pj = zeros(1,n);
    pj(pj<=2*preq) = 1;
    pj = pj/sum(pj);

    phi = (P.*pj)./(P*pj');
    Iold = H(pj) - J(P,pj,phi);
    I = Iold;

    converged = false;
    thresh = 1e-6;
    t = 1;
    while ~converged
        pj = p_up(P,phi,j,preq);
        phi = (P.*pj)./(P*pj');
        I(t) = H(pj) - J(P,pj,phi);
    %     assert(I>=Iold)
        converged = abs(I(t)-Iold) < thresh;
        Iold = I(t);
        t = t+1;
    end

    px(:,k) = pj;
    
    figure(1)
    plot(I)
    figure(2)
    semilogy(j,pj)
    C(k) = I(end)/log(2)*Rs;  
end

%%
figure(3)
plot(PdBm,C/Rs)
grid on

% save(sprintf('apdcap%ddB.mat',GdB),'PdBm','C','power','px','G','R','kA','N0','Rs','df');

function j = J(P,p,phi)
    tmp = P.*log(phi);
    tmp(isinf(tmp))=0;
    tmp(isnan(tmp))=0;
    j = -sum(tmp*p');
end

function p = p_up(P,phi,power,preq)
    tmp = P.*log(phi);
    tmp(isinf(tmp))=0;
    tmp(isnan(tmp))=0;
    r = exp(sum(tmp));   
    
    nu = fzero(@(nu) sum(power.*r.*exp(-nu*power))./sum(r.*exp(-nu*power))-preq,1);
    tmp = r.*exp(-nu*power);
    
    p = tmp/sum(tmp);
end
