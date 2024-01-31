
BPP_ = [];
SSIM_ = [];
ii = 220001;
jj = 220048;
for i = ii:1:jj
    try
    file = load(strcat('PAP7_OMP_SDIC_ALL_file_(',int2str(i),').mat'));
    BPP_   = [BPP_;file.OMP_SDIC_RESULTS_ALL.OVERALL_BPP];
    SSIM_  = [SSIM_;file.OMP_SDIC_RESULTS_ALL.OVERALL_SSIM_FIN];
      
    catch
    end
end
    mod = 30;
    for i = 2:2:30
        ind = find(BPP_(:,1) == i & BPP_(:,2) == mod);
        A = BPP_(ind,:);
        BPP(i/2,:) = mean(BPP_(ind,:),1);
        SSIM(i/2,:) = mean(SSIM_(ind,:),1);
    end
    
    plot(SSIM(:,1) ,BPP(:,4));
    hold on;
    plot(SSIM(:,4), BPP(:,5));
    plot(SSIM(:,7),BPP(:,6));
    plot(SSIM(:,10),BPP(:,7));
    plot(SSIM(:,13),BPP(:,8));
   
    %ylim([0 1]);
    xtitle = 'ssim';
    cc = 0;
    hold off
    xxxx = 0;