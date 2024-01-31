%PAP7_OMP_EOMP_file_(file_index).mat
%        OVERALL_BPP_RESULTS = [OVERALL_BPP_RESULTS; err_threshold(ii), mod_threshold(jj),...
%            BPP_ORI,BPP_OMP,BPP_OMP_ent, BPP_OMP_BRe, BPP_OMP_BRe_ent, ...
%            BPP_EOMP,BPP_EOMP_ent, BPP_EOMP_BRe, BPP_EOMP_BRe_ent];

%        OVERALL_SSIM_RESULTS = [OVERALL_SSIM_RESULTS; err_threshold(ii), mod_threshold(jj),...
%            SSIM_OMP,SSIM_OMP_BRe,SSIM_EOMP,SSIM_EOMP_BRe];

%        OVERALL_PSNR_RESULTS = [OVERALL_PSNR_RESULTS; err_threshold(ii), mod_threshold(jj),...
%            PSNR_OMP,PSNR_OMP_BRe,PSNR_EOMP,PSNR_EOMP_BRe];

%        OVERALL_NMSE_RESULTS = [OVERALL_NMSE_RESULTS; err_threshold(ii), mod_threshold(jj),...
%            NMSE_OMP,NMSE_OMP_BRe,NMSE_EOMP,NMSE_EOMP_BRe];

%PAP7_OTHERS_file_(file_index).mat
%PAP7_OMP_EOMP_file_(file_index).mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SSIM_1 = [];
PSNR_1 = [];
BPP__1 = [];

SSIM_2 = {[],[],[],[],[],[],[],[],[],[],[]};
PSNR_2 = {[],[],[],[],[],[],[],[],[],[],[]};
BPP__2 = {[],[],[],[],[],[],[],[],[],[],[]};
SSIM_2_m = {[],[],[],[],[],[],[],[],[],[],[]};
PSNR_2_m = {[],[],[],[],[],[],[],[],[],[],[]};
BPP__2_m = {[],[],[],[],[],[],[],[],[],[],[]};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
study = 'ssim_psnr_bpp';
study = 'ssim_psnr_bpp_others';

switch study
    case 'ssim_psnr_bpp'
        files_to_be_excluded = [];
        for file_index = 320001:1:322000
            try
                res1 = load(strcat('PAP7_OMP_EOMP_file_(',int2str(file_index),').mat'));
                res2 = load(strcat('PAP7_OTHERS_file_(',int2str(file_index),').mat'));
                SSIM_1 = [SSIM_1;res1.SDIC_RESULTS.OVERALL_SSIM_RESULTS];
                PSNR_1 = [PSNR_1;res1.SDIC_RESULTS.OVERALL_PSNR_RESULTS];
                BPP__1 = [BPP__1;res1.SDIC_RESULTS.OVERALL_BPP_RESULTS];

                SSIM_2_1 = [SSIM_2_1;res2.OTHERS_RESULTS_SA.USQ(:,2)];
                BPP__2_1 = [BPP__2_1;res2.OTHERS_RESULTS_SA.USQ(:,1)];
                SSIM_2_2 = [SSIM_2_2;res2.OTHERS_RESULTS_SA.NUSQ(:,2)];
                BPP__2_2 = [BPP__2_2;res2.OTHERS_RESULTS_SA.NUSQ(:,1)];
                SSIM_2_3 = [SSIM_2_3;res2.OTHERS_RESULTS_SA.VQ(:,2)];
                BPP__2_3 = [BPP__2_3;res2.OTHERS_RESULTS_SA.VQ(:,1)];
                SSIM_2_4 = [SSIM_2_4;res2.OTHERS_RESULTS_SA.gbl_mmc_h(:,2)];
                BPP__2_4 = [BPP__2_4;res2.OTHERS_RESULTS_SA.gbl_mmc_h(:,1)];
                SSIM_2_5 = [SSIM_2_5;res2.OTHERS_RESULTS_SA.gbl_mmc_f(:,2)];
                BPP__2_5 = [BPP__2_5;res2.OTHERS_RESULTS_SA.gbl_mmc_f(:,1)];
                SSIM_2_6 = [SSIM_2_6;res2.OTHERS_RESULTS_SA.lvl_mmc(:,2)];
                BPP__2_6 = [BPP__2_6;res2.OTHERS_RESULTS_SA.lvl_mmc(:,1)];
                SSIM_2_7 = [SSIM_2_7;res2.OTHERS_RESULTS_SA.spiht(:,2)];
                BPP__2_7 = [BPP__2_7;res2.OTHERS_RESULTS_SA.spiht(:,1)];
                SSIM_2_8 = [SSIM_2_8;res2.OTHERS_RESULTS_SA.ezw(:,2)];
                BPP__2_8 = [BPP__2_8;res2.OTHERS_RESULTS_SA.ezw(:,1)];
                SSIM_2_9 = [SSIM_2_9;res2.OTHERS_RESULTS_SA.wdr(:,2)];
                BPP__2_9 = [BPP__2_9;res2.OTHERS_RESULTS_SA.wdr(:,1)];
            catch
                files_to_be_excluded = [files_to_be_excluded; file_index];
            end

            check = 0;


        end
        for i =1:1:50
            SSIM_1_m(i,:) = mean(SSIM_1(i:50:end,:),1);
            PSNR_1_m(i,:) = mean(PSNR_1(i:50:end,:),1);
            BPP__1_m(i,:) = mean(BPP__1(i:50:end,:),1);
            check = 0;
        end

        subplot(1,2,1)
        hold on;
        plot(BPP__1_m(1:5:50,4),SSIM_1_m(1:5:50,3),'LineWidth',2);%OMP
        %plot(SSIM_1_m(1:5:50,3),BPP__1(1:5:50,5));%OMP_ENT
        plot(BPP__1_m(1:5:50,6),SSIM_1_m(1:5:50,3),'LineWidth',2);%OMP_BRe
        %plot(SSIM_1_m(1:5:50,3),BPP__1(1:5:50,7));%OMP_BRe_ENT
        plot(BPP__1_m(1:5:50,8),SSIM_1_m(1:5:50,9),'LineWidth',2);%EOMP
        %plot(SSIM_1_m(1:5:50,9),BPP__1(1:5:50,9));%EOMP_ENT
        plot(BPP__1_m(1:5:50,10),SSIM_1_m(1:5:50,9),'LineWidth',2);%EOMP_BRe
        %plot(SSIM_1_m(1:5:50,9),BPP__1(1:5:50,11));%EOMP_BRe_ENT
        legend({'OMP','OMP-BRe','EOMP','EOMP-BRe'},"Location",'southeast')
        xlabel('BPP');
        ylabel('SSIM');
        title('SDic, \alpha = 10')
        grid on;
        %hold on;
        %plot(SSIM_2_(:,end),BPP__2(:,end));
        set(gca,'FontSize',12,'fontWeight','bold')
        hold off;
        subplot(1,2,2)
        hold on;
        plot(BPP__1_m(1:5:50,4),PSNR_1_m(1:5:50,3),'LineWidth',2);%OMP
        %plot(SSIM_1_m(1:5:50,3),BPP__1(1:5:50,5));%OMP_ENT
        plot(BPP__1_m(1:5:50,6),PSNR_1_m(1:5:50,3),'LineWidth',2);%OMP_BRe
        %plot(SSIM_1_m(1:5:50,3),BPP__1(1:5:50,7));%OMP_BRe_ENT
        plot(BPP__1_m(1:5:50,8),PSNR_1_m(1:5:50,9),'LineWidth',2);%EOMP
        %plot(SSIM_1_m(1:5:50,9),BPP__1(1:5:50,9));%EOMP_ENT
        plot(BPP__1_m(1:5:50,10),PSNR_1_m(1:5:50,9),'LineWidth',2);%EOMP_BRe
        %plot(SSIM_1_m(1:5:50,9),BPP__1(1:5:50,11));%EOMP_BRe_ENT
        legend({'OMP','OMP-BRe','EOMP','EOMP-BRe'},"Location",'southeast');
        xlabel('BPP');
        ylabel('PSNR');
        title('SDic, \alpha = 10')
        grid on;
        %hold on;
        %plot(SSIM_2_(:,end),BPP__2(:,end));
        set(gca,'FontSize',12,'fontWeight','bold')
        check = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'ssim_psnr_bpp_others'
        files_to_be_excluded = [];
        for file_index = 321001:1:322000
            try
                %res1 = load(strcat('PAP7_FBP_EFBP_file_(',int2str(file_index),').mat'));
                res1 = load(strcat('PAP7_OMP_EOMP_file_(',int2str(file_index),').mat'));
                %res1 = load(strcat('PAP7_OOMP_EOOMP_file_(',int2str(file_index),').mat'));
                res2 = load(strcat('PAP7_OTHERS_file_(',int2str(file_index),').mat'));
                SSIM_check = res2.OTHERS_RESULTS_SA.JPG(:,2);
                PSNR_check = res2.OTHERS_RESULTS_SA.JP2(:,3);

                SSIM_1 = [SSIM_1;res1.SDIC_RESULTS.OVERALL_SSIM_RESULTS];
                PSNR_1 = [PSNR_1;res1.SDIC_RESULTS.OVERALL_PSNR_RESULTS];
                BPP__1 = [BPP__1;res1.SDIC_RESULTS.OVERALL_BPP_RESULTS];

                SSIM_2{1} = [SSIM_2{1};res2.OTHERS_RESULTS_SA.JPG(:,2)];
                PSNR_2{1} = [PSNR_2{1};res2.OTHERS_RESULTS_SA.JPG(:,3)];
                BPP__2{1} = [BPP__2{1};res2.OTHERS_RESULTS_SA.JPG(:,1)];

                SSIM_2{2} = [SSIM_2{2};res2.OTHERS_RESULTS_SA.JP2(:,2)];
                PSNR_2{2} = [PSNR_2{2};res2.OTHERS_RESULTS_SA.JP2(:,3)];
                BPP__2{2} = [BPP__2{2};res2.OTHERS_RESULTS_SA.JP2(:,1)];

                SSIM_2{3} = [SSIM_2{3};res2.OTHERS_RESULTS_SA.USQ(:,2)];
                PSNR_2{3} = [PSNR_2{3};res2.OTHERS_RESULTS_SA.USQ(:,3)];
                BPP__2{3} = [BPP__2{3};res2.OTHERS_RESULTS_SA.USQ(:,1)];

                SSIM_2{4} = [SSIM_2{4};res2.OTHERS_RESULTS_SA.NUSQ(:,2)];
                PSNR_2{4} = [PSNR_2{4};res2.OTHERS_RESULTS_SA.NUSQ(:,3)];
                BPP__2{4} = [BPP__2{4};res2.OTHERS_RESULTS_SA.NUSQ(:,1)];

                SSIM_2{5} = [SSIM_2{5};res2.OTHERS_RESULTS_SA.VQ(:,2)];
                PSNR_2{5} = [PSNR_2{5};res2.OTHERS_RESULTS_SA.NUSQ(:,3)];
                BPP__2{5} = [BPP__2{5};res2.OTHERS_RESULTS_SA.VQ(:,1)];

                SSIM_2{6} = [SSIM_2{6};res2.OTHERS_RESULTS_SA.gbl_mmc_h(:,2)];
                PSNR_2{6} = [PSNR_2{6};res2.OTHERS_RESULTS_SA.gbl_mmc_h(:,3)];
                BPP__2{6} = [BPP__2{6};res2.OTHERS_RESULTS_SA.gbl_mmc_h(:,1)];

                SSIM_2{7} = [SSIM_2{7};res2.OTHERS_RESULTS_SA.gbl_mmc_f(:,2)];
                PSNR_2{7} = [PSNR_2{7};res2.OTHERS_RESULTS_SA.gbl_mmc_f(:,3)];
                BPP__2{7} = [BPP__2{7};res2.OTHERS_RESULTS_SA.gbl_mmc_f(:,1)];

                SSIM_2{8} = [SSIM_2{8};res2.OTHERS_RESULTS_SA.lvl_mmc(:,2)];
                PSNR_2{8} = [PSNR_2{8};res2.OTHERS_RESULTS_SA.lvl_mmc(:,3)];                
                BPP__2{8} = [BPP__2{8};res2.OTHERS_RESULTS_SA.lvl_mmc(:,1)];

                SSIM_2{9} = [SSIM_2{9};res2.OTHERS_RESULTS_SA.spiht(:,2)];
                PSNR_2{9} = [PSNR_2{9};res2.OTHERS_RESULTS_SA.spiht(:,3)];
                BPP__2{9} = [BPP__2{9};res2.OTHERS_RESULTS_SA.spiht(:,1)];

                SSIM_2{10} = [SSIM_2{10};res2.OTHERS_RESULTS_SA.ezw(:,2)];
                PSNR_2{10} = [PSNR_2{10};res2.OTHERS_RESULTS_SA.ezw(:,3)];
                BPP__2{10} = [BPP__2{10};res2.OTHERS_RESULTS_SA.ezw(:,1)];

                SSIM_2{11} = [SSIM_2{11};res2.OTHERS_RESULTS_SA.wdr(:,2)];
                PSNR_2{11} = [PSNR_2{11};res2.OTHERS_RESULTS_SA.wdr(:,3)];
                BPP__2{11} = [BPP__2{11};res2.OTHERS_RESULTS_SA.wdr(:,1)];
                
            catch
                files_to_be_excluded = [files_to_be_excluded; file_index];
            end

            check = 0;

        end
        for i =1:1:150
            SSIM_1_m(i,:) = mean(SSIM_1(i:150:end,:),1);
            PSNR_1_m(i,:) = mean(PSNR_1(i:150:end,:),1);
            BPP__1_m(i,:) = mean(BPP__1(i:150:end,:),1);
            check = 0;
        end
        bpp_int = [0:0.2:3];
        for i = 1:1:15
            for j = 1:1:11
                ind = find (BPP__2{j}<bpp_int(i+1) & BPP__2{j}>bpp_int(i));
                BPP__2_m{j} = [BPP__2_m{j};mean(BPP__2{j}(ind))];
                SSIM_2_m{j} = [SSIM_2_m{j};mean(SSIM_2{j}(ind))];
                PSNR_2_m{j} = [PSNR_2_m{j};mean(PSNR_2{j}(ind))];
                cc = 0;
            end
        end
        ALGS_GROUP_1 = {'OMP','Ent(OMP)','OMP-BRe','Ent(OMP-BRe)','EOMP','Ent(EOMP)','EOMP-BRe','Ent(EOMP-BRe)'};
        ALGS_GROUP_2 = {'JPG','JP2','USQ','NUSQ','VQ','gbl_mmc_h','gbl_mmc_f','lvl_mmc','spiht','ezw','wdr'};
        AG1_ssim_ind = {[3,4,5],[3,4,5],[6,7,8],[6,7,8],[9,10,11],[9,10,11],[12,13,14],[12,13,14]};
        AG1_psnr_ind = {[3,4,5],[3,4,5],[6,7,8],[6,7,8],[9,10,11],[9,10,11],[12,13,14],[12,13,14]};
        AG1_nmse_ind = {[3,4,5],[3,4,5],[6,7,8],[6,7,8],[9,10,11],[9,10,11],[12,13,14],[12,13,14]};
        AG1_bpp__ind = {[4],[5],[6],[7],[8],[9],[10],[11]};
       % OMP           1,     JPG       1
       % Ent(OMP)      2,     JP2       2
       % OMP-BRe       3,     USQ       3
       % Ent(OMP-BRe)  4,     NUSQ      4
       % EOMP          5,     VQ        5
       % Ent(EOMP)     6,     gbl_mmc_h 6
       % EOMP-BRe      7,     gbl_mmc_f 7
       % Ent(EOMP-BRe) 8,     lvl_mmc   8
       %                      spiht     9            
       %                      ezw       10
       %                      wdr       11
       %
        sel_mod = [1,2,3,4,5];%10,20,30,40,50
        sel_m = 1;
        sel_AG1 = [1,2,4,6,8];
        sel_AG2 = [9,10,11];
        subplot(1,2,1)
        hold on;
        for i = 1:1:size(sel_AG1,2) 
            plot(BPP__1_m(sel_mod(sel_m):15:150,AG1_bpp__ind{sel_AG1(i)}),...
                SSIM_1_m(sel_mod(sel_m):15:150,AG1_ssim_ind{sel_AG1(i)}(2))...%(1) refers to no-denoise
                ,'LineWidth',2);
        end
        
        for i = 1:1:size(sel_AG2,2)
            plot(BPP__2_m{sel_AG2(i)},SSIM_2_m{sel_AG2(i)},'LineWidth',2);
            %plot(0:0.01:3,polyval(polyfit(BPP__2{i},SSIM_2{i},4),0:0.01:3),'LineWidth',2);      
        end     
        ALGS = {ALGS_GROUP_1{sel_AG1}, ALGS_GROUP_2{sel_AG2}};
        legend(ALGS,"Location",'southeast');
        xlabel('BPP');   ylabel('SSIM');    title('SDic, \alpha = 10');
        grid on;
        set(gca,'FontSize',12,'fontWeight','bold');
        hold off;
        
        
        subplot(1,2,2)
        hold on;
        for i = 1:1:size(sel_AG1,2) 
            plot(BPP__1_m(sel_mod(sel_m):15:150,AG1_bpp__ind{sel_AG1(i)}),...
                PSNR_1_m(sel_mod(sel_m):15:150,AG1_psnr_ind{sel_AG1(i)}(2))...%(1) refers to no-denoise
                ,'LineWidth',2);
        end
        
        for i = 1:1:size(sel_AG2,2)
            plot(BPP__2_m{sel_AG2(i)},PSNR_2_m{sel_AG2(i)},'LineWidth',2);
            %plot(0:0.01:3,polyval(polyfit(BPP__2{i},SSIM_2{i},4),0:0.01:3),'LineWidth',2);      
        end     
        ALGS = {ALGS_GROUP_1{sel_AG1}, ALGS_GROUP_2{sel_AG2}};
        legend(ALGS,"Location",'southeast');
        xlabel('BPP');   ylabel('PSNR');    title('SDic, \alpha = 10');
        grid on;
        set(gca,'FontSize',12,'fontWeight','bold');
        hold off;
        check = 0;
        %imbinarize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end