% APDBW="Inf";
for APDBW = ["Inf", "25", "10"]
    for ModBW = ["25" "10"]
        for p = 2:3
            for L= ["2" "7"]
                pam_access_nets_qsub(num2str(2^p),L,"1550",ModBW,"apd",APDBW,"15","5","2");
            end
        end
    end
end