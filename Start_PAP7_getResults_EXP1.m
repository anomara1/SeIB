img1_res = load("PAP7_OMP_EOMP_file_(321004).mat");
%img2_res = load("PAP7_OMP_EOMP_file_(220046).mat");

subplot(1,3,1)
for i = 10:10:50
    for alpha = 1:1:15
        ALL1 = img1_res.SDIC_RESULTS.ALL_STHRE_MTHRE_NSEBT;
        ind = find(ALL1(:,1) == i & ALL1(:,2) ==  alpha *2);
        k_alg(alpha,1) = mean(ALL1(ind,3));
        k_enh_alg(alpha,1) = mean(ALL1(ind,10));
        
        
    end
    %to_plot = [[10:10:50]',k_alg,k_enh_alg];
    plot([2:2:30],k_enh_alg./k_alg);
    hold on;
    %plot([10:10:50],k_enh_alg);
end
legend({'Err_{th} = 10%','Err_{th} = 20%','Err_{th} = 30%','Err_{th} = 40%','Err_{th} = 50%'})
xlabel('\alpha'); ylabel('k_{EOMP}/k_{OMP}')
cc = 0;
hold off;
subplot(1,3,2)
for i = 10:10:50
    for alpha = 1:1:15
        ALL1 = img1_res.SDIC_RESULTS.ALL_STHRE_MTHRE_NSEBT;
        ind = find(ALL1(:,1) == i & ALL1(:,2) ==  alpha *2);
        t_alg(alpha,1) = mean(ALL1(ind,9));
        t_enh_alg(alpha,1) = mean(ALL1(ind,16));
        
        
    end
    %to_plot = [[10:10:50]',t_alg,t_enh_alg];
    plot([2:2:30],t_enh_alg./t_alg);
    hold on;
    %plot([10:10:50],k_enh_alg);
end
legend({'Err_{th} = 10%','Err_{th} = 20%','Err_{th} = 30%','Err_{th} = 40%','Err_{th} = 50%'})
xlabel('\alpha'); ylabel('t_{EOMP}/t_{OMP}')
cc = 0;

hold off;
subplot(1,3,3)

for i = 10:10:50
    for alpha = 1:1:15
        ALL1 = img1_res.SDIC_RESULTS.ALL_STHRE_MTHRE_NSEBT;
        ind = find(ALL1(:,1) == i & ALL1(:,2) ==  alpha *2);
        k_alg = ALL1(ind,3);
        k_enh_alg = ALL1(ind,10);
        pro(alpha,1) = sum(k_enh_alg>k_alg)/(size(ind,1));
        
    end
    %to_plot = [[10:10:50]',k_alg,k_enh_alg];
    plot([2:2:30],pro);
    hold on;
    %plot([10:10:50],k_enh_alg);
end
legend({'Err_{th} = 10%','Err_{th} = 20%','Err_{th} = 30%','Err_{th} = 40%','Err_{th} = 50%'})
xlabel('\alpha'); ylabel('P_{FE}')
cc = 0;