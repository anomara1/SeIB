bpp1 = [];
bpp2 = [];
bpp3 = [];
ss1 = [];
ss2 = [];
ss3 = [];
ref_b = [];
ref_s = [];
ii = 321001;
jj = 321150;
for i = ii:1:jj
    file = load(strcat('PAP7_OMP_EOMP_ACC_IAM_IBM_file_(',int2str(i),').mat'));
    %v2 = load(strcat('PAP7_OMP_EOMP_VAR2_file_(',int2str(i),').mat'));
    %v3 = load(strcat('PAP7_OMP_EOMP_VAR3_file_(',int2str(i),').mat'));
    ref = file.SDIC_RESULTS_ACC.OVERALL_BPP_RESULTS(:,[1 2]);

    ref_b = [ref_b , file.SDIC_RESULTS_ACC.OVERALL_BPP_RESULTS(:,4)];
    ref_s = [ref_s, file.SDIC_RESULTS_ACC.OVERALL_SSIM_RESULTS(:,3)];

    bpp1 = [bpp1 , file.SDIC_RESULTS_ACC.OVERALL_BPP_RESULTS(:,8)];
    bpp2 = [bpp2 , file.SDIC_RESULTS_IAM.OVERALL_BPP_RESULTS(:,8)];
    bpp3 = [bpp3 , file.SDIC_RESULTS_IBM.OVERALL_BPP_RESULTS(:,8)];
    ss1 = [ss1 , file.SDIC_RESULTS_ACC.OVERALL_SSIM_RESULTS(:,9)];
    ss2 = [ss2 , file.SDIC_RESULTS_IAM.OVERALL_SSIM_RESULTS(:,9)];
    ss3 = [ss3 , file.SDIC_RESULTS_IBM.OVERALL_SSIM_RESULTS(:,9)];
    

end
    ref_b_ = nanmean(ref_b,2);
    ref_s_ = nanmean(ref_s,2);
    bpp1_ = nanmean(bpp1,2);
    bpp2_ = nanmean(bpp2,2);
    bpp3_ = nanmean(bpp3,2);
    ss1_ = nanmean(ss1,2);
    ss2_ = nanmean(ss2,2);
    ss3_ = nanmean(ss3,2);
    ind = find(ref(:,1) == 10);
    
    %plot(ref_s(ind,1)./(jj-ii+1),ref_b_(ind,1)./ref_b_(ind,1) );
    plot(ref_s_(ind,1),ref_b_(ind,1));
    hold on;

    plot(ss1_(ind,1),bpp1_(ind,1));%./ref_b_(ind,1) );
    plot(ss2_(ind,1),bpp2_(ind,1));%./ref_b_(ind,1));
    plot(ss3_(ind,1),bpp3_(ind,1));%./ref_b_(ind,1));
    %ylim([0 1]);
    xtitle = 'ssim';
    cc = 0;
    hold off
    xxxx = 0;