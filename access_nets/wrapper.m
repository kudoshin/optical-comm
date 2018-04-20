% % APDBW="Inf";
% for APDBW = [Inf, 25, 10]
%     for ModBW = [25, 10]
%         for p = 2:3
%             for L= [2 7]
%                 pam_access_nets_qsub(num2str(2^p),L,"1550",ModBW,"apd",APDBW,"15","5","2");
%             end
%         end
%     end
% end

%%
L = [40 50];

ModBW = 25;
APDBW = 25;
M = [2,4,8];
for l=L
    for m=M
        pam_access_nets_qsub(m, l, 1550, ModBW, 0, Inf, 15, 5, 1);
    end
end

%%
% pam_access_nets_qsub(8,0,1550,25,2,25,15,5,1);
