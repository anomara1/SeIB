function  runAVGs()
RES_BBB_BIG_DS_128 = load('BBB_BIG_DATASET_RESULTS');
GI_1 = [2 3 4 5 12 13 14 16 21 25 27 29 31 36 40 41 46 52  55  56 ... 
    57 60  61 62 64  67 74 75 76 77  78   80 81  82 83 84  85 88  94 95 90 93 ...
    97 98  99  100];
RES_BBB_DS_128 = load('BBB_DATASET_RESULTS_128_128');
%(1039,1033,1031,1027,1017,1016,1015,1014,1013,1011,1010,1007,1006,1001 semi)
GI_2 = [38     39 33 31 27 17 16 15 14 13 11 10 7 6];
RES_BBB_DS_256 = load('BBB_DATASET_RESULTS_256_256');
GI_3 = [38     39 33 31 27 17 16 15 14 13 11 10 7 6];
RES_BBB_DS_512 = load('BBB_DATASET_RESULTS_512_512');
GI_4 = [38     39 33 31 27 17 16 15 14 13 11 10 7 6];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RES_MiOMP_AAA_DS_256 = load('AAA_MiOMP_40_256_256');
GI_5 = [14 16 18 24 27 29 32 35 37 40 41 42 43 44];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res_ =RES_MiOMP_AAA_DS_256; %RES_BBB_DS_128;%RES_BBB_DS_256;
GI_sel = GI_5;%1:40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = 0;
res_ssim_fixed = [];
res_ssim_fixed_temp = [];
res_psnr_fixed = [];
res_psnr_fixed_temp = [];
res_ttest = [];
res_ttest_temp = [];
res_backward =[];
res_backward_temp =[];
res_forward =[];
res_forward_temp =[];

for i = 1:1:10
    for j =1:1:length(GI_sel)%100
        GI = GI_sel;
       count = 1;
       X = res_.RES(j);%{GI(j)};%.res;
       res_ssim_fixed = [res_ssim_fixed; X.SSIM_Fixed_sp_res(i,:)];
       res_psnr_fixed = [res_psnr_fixed ; X.PSNR_FIXED_sp_res(i,:)];
       res_ttest = X.SSIM_ttest;
        res_backward = [res_backward ; X.Backward_Algs(i,:)];
        res_forward = [res_forward ; X.Forward_Algs(i,:)];       
       cc=0;
    end
    ss = res_ssim_fixed;
    
    res_ssim_fixed_temp(i,:)=nanmean(res_ssim_fixed);
    res_ttest(i,:) = [ res_ssim_fixed_temp(i,1) my_ttest(ss(:,2),ss(:,3)) my_ttest(ss(:,2),ss(:,4)) my_ttest(ss(:,2),ss(:,5)) my_ttest(ss(:,2),ss(:,6)) ...
        my_ttest(ss(:,2),ss(:,7)) my_ttest(ss(:,2),ss(:,8)) my_ttest(ss(:,2),ss(:,9)) my_ttest(ss(:,2),ss(:,10)) my_ttest(ss(:,2),ss(:,11))];
    res_psnr_fixed_temp(i,:)=nanmean(res_psnr_fixed);
    res_backward_temp(i,:)=nanmean(res_backward);
    res_forward_temp(i,:)=nanmean(res_forward);
    res_ssim_fixed = [];
    res_psnr_fixed = [];
    res_backward = [];
    res_forward = [];
    cc = 0;
   

   % cc = 0;
end
cc = 0;
end
