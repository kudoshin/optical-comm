
% assignments = perms(1:7);
% 
% minp = 1;
% 
% for x=assignments'
%     pb = err_pr([0; x],0.2);
%     if pb < minp
%         minp = pb;
%         opt_assign = [0; x];
%     end
% end

%%
% pb = err_pr([7 6 5 1 0 4 2 3],0.2);  % pb=0.0233
% pb = err_pr([0 2 6 1 3 7 4 5],0.2);  % pb=0.0186

sig_const = [0 0;0 1;0 2;1 0;1 1;1 2;2 0;2 1];
sigma = fzero(@(x) 0.01 - err_pr([0 2 6 1 3 7 4 5],x),1/qfuncinv(0.01/2)/2);
SNR = 10*log10(1/sigma^2*sum(sig_const(:).^2)); % SNR = 28.1 dB
sig_const = [0 0;0 1;0 2;1 0;1 2;2 0;2 1; 2 2];
sigma = fzero(@(x) 0.01 - err_pro([0 1 3 4 2 5 7 6],x),1/qfuncinv(0.01/2)/2);
SNRsq = 10*log10(1/sigma^2*sum(sig_const(:).^2)); % SNR = 28.7 dB

%% I(X;Y)
sigma=linspace(0,1);
IXY = zeros(size(sigma));
for i=1:length(sigma)
    Q1 = qfunc(1/2/sigma(i));
    Q2 = qfunc(3/2/sigma(i));
    ptri = 0;
    pxy = [(1-Q1)^2 (Q1-Q2)*(1-Q1) Q2*(1-Q1) (Q1-Q2)*(1-Q1) (Q1-Q2)^2 Q2*(Q1-Q2/2) Q2*(1-Q1) Q2*(Q1-Q2/2)
        Q1*(1-Q1) (1-2*Q1)*(1-Q1) Q1*(1-Q1) Q1*(Q1-Q2) (Q1-Q2)*(1-2*Q1) Q1*(Q1-Q2/2) Q2*Q1 Q2*(1-3/2*Q1)
        Q2*(1-Q1) (Q1-Q2)*(1-Q1) (1-Q1)^2 Q2*(Q1-Q2) (Q1-Q2)^2 (Q1-Q2/2)*(1-Q1) Q2^2 (Q1-Q2+(1-Q1)/2)*Q2
        Q1*(1-Q1) Q1*(Q1-Q2) Q2*Q1 (1-Q1)*(1-2*Q1) (1-2*Q1)*(Q1-Q2) (1-3/2*Q1)*Q2 Q1*(1-Q1) (Q1-Q2/2)*Q1
        Q1^2 Q1*(1-2*Q1) Q1^2 Q1*(1-2*Q1) (1-2*Q1)^2 Q1*(1-3/2*Q1) Q1^2 Q1*(1-3/2*Q1)
        Q2*Q1 Q1*(Q1-Q2) Q1*(1-Q1) Q2*(1-2*Q1) (Q1-Q2)*(1-2*Q1) (1-Q1)*(1-Q1)-ptri Q2*Q1 (Q1-Q2)*Q1+ptri
        Q2*(1-Q1) (Q1-Q2)*Q2 Q2^2 (Q1-Q2)*(1-Q1) (Q2-Q1)^2 (Q1-Q2+(1-Q1)/2)*Q2 (1-Q1)^2 (1-Q1)*(Q1-Q2/2)
        Q1*Q2 Q2*(1-2*Q1) Q1*Q2 (Q1-Q2)*Q1 (Q1-Q2)*(1-2*Q1) (Q1-Q2)*Q1+ptri Q1*(1-Q1) (1-Q1)*(1-Q1)-ptri]; 

    py = mean(pxy);
    HY = -sum(py.*log2(py));
    HYX = -mean(sum(pxy.*log2(pxy),2));
    IXY(i) = HY-HYX;
end
SNRdB = 10*log10(9/8/sigma.^2);
plot(SNRdB,IXY/2)
%%
% s2
% ^
% | x x
% | x x x     2 symbols - 3 bits
% | x x x
% -------> s1

function pb = err_pr(bit_assign,sigma)
d = 3*squareform(pdist(de2bi(bit_assign),'hamming'));
Q = qfunc(1/2/sigma);
x = -1/2:0.001:3;
p86 = Q^2 + 1/sqrt(2*pi*sigma^2)*trapz(x,exp(-x.^2/(2*sigma^2)).*qfunc((x+1)/sigma));

pb  = 1/4*((Q-Q^2)*(d(1,2)+d(1,4)+d(7,4)+d(2,3)+d(3,6)+d(7,8))...
    + Q^2*(d(1,5)+d(2,4)+d(2,6)+d(3,5)+d(4,8)+d(5,7)+d(5,8)/2+d(5,6)/2)...
    + (Q-2*Q^2)*(d(2,5)+d(4,5)+d(5,6)+d(5,8))...
    + p86*d(6,8));
% pes = [(Q-Q^2)*(d(1,2)+d(1,4))+Q^2*d(1,5)
%     (Q-Q^2)*(d(1,2)+d(2,3))+Q^2*(d(3,4)+d(3,6))+(Q-2*Q^2)*d(2,5)
%     (Q-Q^2)*(d(2,3)+d(3,6))+Q^2*d(3,5)
%     (Q-Q^2)*(d(1,4)+d(4,7))+Q^2*(d(4,8)+d(2,4))+(Q-2*Q^2)*d(4,5)
%     (Q-2*Q^2)*(d(2,5)+d(4,5)+d(5,6)+d(5,8))+Q^2*(d(1,5)+d(3,5)+d(5,7)+(d(5,6)+d(5,8))/2)
%     (Q-Q^2)*d(3,6)+Q^2*d(2,6)+(Q-2*Q^2)*d(5,6)+p86*d(6,8)
%     (Q-Q^2)*(d(4,7)+d(7,8))+Q^2*d(5,7)
%     (Q-Q^2)*(d(5,8)+d(7,8))+Q^2*d(4,8)+p86*d(6,8)];

% sig_const = [0 0;0 1;0 2;1 0;1 1;1 2;2 0;2 1];
% distances = pdist(sig_const);
% NNi = distances == 1;
% NNdiagi = distances == sqrt(2);
% d1 = hist(d(NNi));
% d2 = hist(d(NNdiagi));
end

%%
% s2
% ^
% | x x x
% | x   x     2 symbols - 3 bits
% | x x x
% -------> s1

function pb = err_prsq(bit_assign, sigma)
d = 3*squareform(pdist(de2bi(bit_assign),'hamming'));
Q = qfunc(1/2/sigma);
x = 0.5:0.0001:1;
pmc = Q/2 - 1/sqrt(2*pi*sigma^2)*trapz(x,exp(-x.^2/(2*sigma^2)).*(1/2-qfunc((1-x)/sigma)));
pes = [(Q - Q^2/2)*(d(1,2)+d(1,4))
    (Q-Q^2)*(d(1,2)+d(2,3))+pmc*(d(2,4)+d(2,5))
    (Q-Q^2/2)*(d(2,3)+d(3,5))
    (Q-Q^2)*(d(1,4)+d(4,6))+pmc*(d(2,4)+d(4,7))
    (Q-Q^2)*(d(3,5)+d(5,8))+pmc*(d(2,5)+d(5,7))
    (Q-Q^2/2)*(d(4,6)+d(6,7))
    (Q-Q^2)*(d(6,7)+d(7,8))+pmc*(d(4,7)+d(5,7))
    (Q-Q^2/2)*(d(5,8)+d(7,8))];
pb = mean(pes);
end
