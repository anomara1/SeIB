%imDim,imName,imExt,dicType,ssimDIF, enh
function Start_PAP7_SISM_v_4(file_index)
clc;
batchSize = 256;%16,8
sparse = 128;
imDim = 32; %512,256,128
imExt ='.jpg'; %jpg, gif,jpeg
try
for i = file_index:file_index%322000
    %t1 = clock;
    ff = load('images_sizes.mat');
    if file_index>320000
        ori_bpp = ff.size_(file_index-320000);
    else
        ori_bpp = 0;
    end
    fName = strcat('(',int2str(i),')',imExt);
    [bytes, jpg_ssim, jpg_psnr, qual] = getJPEG_qualme_ano(fName);
    JPEG_RES = [bytes*8/(256*256), jpg_ssim, jpg_psnr, qual];
    [bytes, j2000_ssim, j2000_psnr, c_rat] = getJPEG2000_qual_me_ano(fName);
    JPEG2000_RES = [bytes*8/(256*256), j2000_ssim, j2000_psnr, c_rat];
    pre_res  = strcat('PAP7_ALL_VARS_RESULTS/PAP7_OMP_SDIC_ALL');
    dicType_ = 'SDic';
    [OMP_SDIC_RESULTS_ALL] = PAP7_main_code(fName,ori_bpp,batchSize,sparse,dicType_,JPEG_RES,JPEG2000_RES,'ALL');
    %[SDIC_RESULTS_IAM] = PAP7_main_code(fName,ori_bpp,batchSize,sparse,dicType_,JPEG_RES,JPEG2000_RES,'IAM');
    %[SDIC_RESULTS_IBM] = PAP7_main_code(fName,ori_bpp,batchSize,sparse,dicType_,JPEG_RES,JPEG2000_RES,'IBM');
    
    dicType_ = 'LDic';
    fName = strcat('(',int2str(i),')',imExt);
    %[LEARNED_RESULTS_SA] = AAA_Second_BE_Final_SparseALGs_enh_PAP5(fName,batchSize,sparse,dicType_);
    LDIC_RESULTS_VAR2 = [];
    %[LDIC_RESULTS_VAR2] = PAP7_main_code(fName,ori_bpp,batchSize,sparse,dicType_,JPEG_RES,JPEG2000_RES);
    
    
    save(strcat(pre_res,'_file_','(',int2str(i),')'),'OMP_SDIC_RESULTS_ALL');
    check = 0;
    
end
catch 
    fid = fopen( strcat('Avoid_(',int2str(file_index),')_.txt'), 'wt' );
    fclose(fid);
end
end
function [RESULTS] = PAP7_main_code(img,ori_bpp,batchSize,sparse_,dicType,JPEG_RES,JPEG2000_RES,IAM_IBM_NOTIFY) %#ok<*INUSD>
%BRE_IMG Summary of this function goes here
%   Detailed explanation goes here
global bits_per_coeff;
global bits_per_index;
global bits_per_NoCoeff;
bits_per_coeff = 9;
bits_per_index = 9;
bits_per_NoCoeff = 7;
global sparse;
sparse = sparse_;
img_pixels = imread(img);
vecLength = batchSize;%batchSize*batchSize;
img_double = im2double(img_pixels);
Dic = getDic_(dicType,vecLength);
[VecMatrix,kr,kc] = AAA_getIMGvectors_(img_double,batchSize);
index = 1;
warning('off','all');

err_threshold = 2:2:30;
mod_threshold = 2:2:30;
OVERALL_TIME = [];
OVERALL_PARA = [];
OVERALL_SSIM_avg = [];
OVERALL_SSIM_FIN = [];
OVERALL_NMSE_avg = [];
OVERALL_NMSE_FIN = [];
OVERALL_PSNR_FIN = [];

OVERALL_BPP  = [];
OVERALL_ALPHAS = [];
TIME_RES = [];
PARA_RES = [];
SSIM_RES = [];
NMSE_RES = [];
BPP_RES  = [];
Alphas = [];

indicator_ = 1;
DicCoh = Dic'*Dic;
%omp_ci_matrix = zeros(size(Dic,2),kr*kc,length(err_threshold));
%eomp_ci_matrix = zeros(size(Dic,2),kr*kc,length(err_threshold),length(mod_threshold));
for ii = 1:1:length(err_threshold)
    do_one_time = 0;
    kk = 1;
    
    
    for jj = 1:1:length(mod_threshold)
        notify = 0;
        if do_one_time == 0
            X_OMP = [];
            OMP_COEFF = [];
            OMP_IND = [];
            Xhat_OMP = [];
            N1 = zeros(kr,kc,length(err_threshold));
            S1  = zeros(kr,kc,length(err_threshold));
            E1  = zeros(kr,kc,length(err_threshold));
            BPP1 = zeros(kr,kc,length(err_threshold));
            T1 = zeros(kr,kc,length(err_threshold));
        end
        PKbM_COEFF = [];
        PKbM_IND = [];
        PPKbM_1_COEFF = [];
        PPKbM_1_IND = [];
        PPKbM_2_COEFF = [];
        PPKbM_2_IND = [];
        UDbM_COEFF = [];
        UDbM_IND = [];
        
        clc;
        fprintf(strcat(img,' ... ',dicType,'_',IAM_IBM_NOTIFY, ' ... ', '(Ours Alg. Reached %.2f%%)\n'),(indicator_-1)*100/(length(err_threshold)*length(mod_threshold)));%2012
        indicator_ = indicator_ + 1;
        
        for i=1:1:kr
            k = 1;
            for j=1:1:kc
                Vector =  im2double(img_pixels(256*(i-1)+1:i*256,j)); %VecMatrix(i,j,:);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(OMP(1990) - StOMP(2006)? - SWOMP(2009)  LASSO(2009) - COSAMP(2008) - GSM(2020)-
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OGA() - OMPvar(2015) - SAEN(2018) - GOMP(2010) -  - cpwwen(2018) - clarswlasso(2018) )%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if do_one_time == 0
                    for v =1:1:kc
                        Vector_ =  im2double(img_pixels(256*(i-1)+1:i*256,v));
                        tic;
                        [~, omp_coeff, omp_iopt, omp_error,omp_ssim, omp_hat]       = AAA_SSIM_OMP_new_pap7_(Dic, Vector_,sparse,1, err_threshold(ii));
                        T1(i,v,ii) = toc;
                        OMP_COEFF = [OMP_COEFF;omp_coeff];
                        OMP_IND = [OMP_IND;omp_iopt];
                        N1(i,v,ii) = size(omp_iopt,1);
                        S1(i,v,ii) = omp_ssim(end,1);
                        E1(i,v,ii) = omp_error(end,1);
                        BPP1(i,v,ii) = (ceil(log2(size(omp_coeff,1))) + bits_per_coeff + (N1(i,v,ii)-1)*(bits_per_coeff + bits_per_index))/256;
                        X_OMP = [X_OMP Vector_];
                        Xhat_OMP = [Xhat_OMP omp_hat(:,end)];
                        
                        kk = kk + 1;
                    end
                    X_PKbM = zeros(size(X_OMP,1),size(X_OMP,2));
                    X_PPKbM_1 = zeros(size(X_OMP,1),size(X_OMP,2));
                    X_PPKbM_2 = zeros(size(X_OMP,1),size(X_OMP,2));
                    X_UDbM = zeros(size(X_OMP,1),size(X_OMP,2));
                end
                if notify == 0
                    if mod(length(Vector),mod_threshold(jj)+1) == 0
                        J_INDEX = [1:mod_threshold(jj):length(Vector)];
                    else
                        J_INDEX = [1:mod_threshold(jj):length(Vector), length(Vector)];
                    end
                    X_PKbM(:,J_INDEX) = Xhat_OMP(:,J_INDEX);
                    X_PPKbM_1(:,J_INDEX) = Xhat_OMP(:,J_INDEX);
                    X_PPKbM_2(:,J_INDEX) = Xhat_OMP(:,J_INDEX);
                    X_UDbM(:,J_INDEX) = Xhat_OMP(:,J_INDEX);
                end
                notify = 1;
                do_one_time = 1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%   Prior Knowledge based Modeling
                %%%%%%%%%%%%%%%%%%%%%%%%   (PKbM)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tic;
                if ismember(j,J_INDEX) == 1
                    PKbM_X_init = zeros(size(Vector,1),1);
                else
                    PKbM_X_init = X_PKbM(:,j-1);
                end
                
                [~, eomp_coeff, eomp_iopt, eomp_error,eomp_SSIM, eomp_hat]       = Enh_Convergence_OMP_pap7_(Dic, Vector,sparse,PKbM_X_init,1, err_threshold(ii));
                PKbM_COEFF = [PKbM_COEFF;eomp_coeff];
                PKbM_IND = [PKbM_IND;eomp_iopt];
                X_PKbM(:,j) = eomp_hat(:,end);
                T2 = toc;
                N2 = size(eomp_iopt,1);
                S2 = eomp_SSIM(end,1);
                E2 = eomp_error(end,1);
                BPP2 = (ceil(log2(size(omp_coeff,1))) + bits_per_coeff + (N2-1)*(bits_per_coeff + bits_per_index))/256;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%   Prior-Posterior Knowledge based
                %%%%%%%%%%%%%%%%%%%%%%%%   Modeling with one parameter
                %%%%%%%%%%%%%%%%%%%%%%%%   (PPKbM_1)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tic;
                if ismember(j,J_INDEX) == 1
                    X_PPKbM_1_init = zeros(size(Vector,1),1);
                    alpha_1 = 0;
                    alpha_2 = 0;
                else
                    [alpha_1, alpha_2, init] = getAlpha(j,J_INDEX,Vector,X_PPKbM_1);
                    X_PPKbM_1_init = init;
                end
                
                [~, eomp_coeff, eomp_iopt, eomp_error,eomp_SSIM, eomp_hat]       = Enh_Convergence_OMP_pap7_(Dic, Vector,sparse,X_PPKbM_1_init,1, err_threshold(ii));
                PPKbM_1_COEFF = [PPKbM_1_COEFF;eomp_coeff];
                PPKbM_1_IND = [PPKbM_1_IND;eomp_iopt];
                X_PPKbM_1(:,j) = eomp_hat(:,end);
                T3 = toc;
                N3 = size(eomp_iopt,1)+2;
                S3 = eomp_SSIM(end,1);
                E3 = eomp_error(end,1);
                BPP3 = (ceil(log2(size(omp_coeff,1))) + bits_per_coeff + (N3-1)*(bits_per_coeff + bits_per_index) - 2* bits_per_index)/256;
                alpha_1_3 = alpha_1;
                alpha_2_3 = alpha_2;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%   Prior-Posterior Knowledge based
                %%%%%%%%%%%%%%%%%%%%%%%%   Modeling with two parameters
                %%%%%%%%%%%%%%%%%%%%%%%%   (PPKbM_2)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tic;
                if ismember(j,J_INDEX) == 1
                    X_PPKbM_2_init = zeros(size(Vector,1),1);
                    alpha_vector(1) = 0;
                    alpha_vector(2) = 0;
                else
                    [alpha_vector, init] = prep_init_matrix(j,J_INDEX,Vector,X_PPKbM_2);
                    X_PPKbM_2_init = init;
                end
                
                [~, eomp_coeff, eomp_iopt, eomp_error,eomp_SSIM, eomp_hat]       = Enh_Convergence_OMP_pap7_(Dic, Vector,sparse,X_PPKbM_2_init,1, err_threshold(ii));
                PPKbM_2_COEFF = [PPKbM_2_COEFF;eomp_coeff];
                PPKbM_2_IND = [PPKbM_2_IND;eomp_iopt];
                X_PPKbM_2(:,j) = eomp_hat(:,end);
                T4 = toc;
                N4 = size(eomp_iopt,1)+2;
                S4 = eomp_SSIM(end,1);
                E4 = eomp_error(end,1);
                BPP4 = (ceil(log2(size(omp_coeff,1))) + bits_per_coeff + (N4-1)*(bits_per_coeff + bits_per_index)-2*bits_per_index)/256;
                alpha_1_4 = alpha_vector(1);
                alpha_2_4 = alpha_vector(2);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%   Prior-Posterior Knowledge based
                %%%%%%%%%%%%%%%%%%%%%%%%   Modeling with Updated Dic
                %%%%%%%%%%%%%%%%%%%%%%%%   (UDbM)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tic;
                X_UDbM_init = zeros(size(Vector,1),1);
                if ismember(j,J_INDEX) == 1
                    Dic_ = Dic;
                    
                else
                    %[Atom_1, Atom_2] = prep_init_matrix(j,J_INDEX,Vector,X_PPKbM_2);
                    [Atom_1,Atom_2] = get_2_bases(j,J_INDEX,Vector,X_UDbM);
                    Dic_ = [Dic,Atom_1, Atom_2];
                end
                
                [~, eomp_coeff, eomp_iopt, eomp_error,eomp_SSIM, eomp_hat]       = Enh_Convergence_OMP_pap7_(Dic_, Vector,sparse,X_UDbM_init,1, err_threshold(ii));
                UDbM_COEFF = [UDbM_COEFF;eomp_coeff];
                UDbM_IND = [UDbM_IND;eomp_iopt];
                X_UDbM(:,j) = eomp_hat(:,end);
                T5 = toc;
                N5 = size(eomp_iopt,1);
                S5 = eomp_SSIM(end,1);
                E5 = eomp_error(end,1);
                BPP5 = (ceil(log2(size(omp_coeff,1))) + bits_per_coeff + (N5-1)*(bits_per_coeff + bits_per_index))/256;
                if size(Dic_,2) > size(Dic,2)
                    
                    alpha_1_5 = eomp_coeff(find(eomp_iopt == size(Dic_,2)-1));
                    alpha_2_5 = eomp_coeff(find(eomp_iopt == size(Dic_,2)));
                    if isempty(alpha_1_5)
                        alpha_1_5 = 0;
                    end
                    
                    
                    if isempty(alpha_2_5)
                        alpha_2_5 = 0;
                    end
                else
                    alpha_1_5 = 0;
                    alpha_2_5 = 0;
                end
                TIME_RES = [TIME_RES;err_threshold(ii), mod_threshold(jj) T1(i,j,ii) T2 T3 T4 T5];
                PARA_RES = [PARA_RES;err_threshold(ii), mod_threshold(jj) N1(i,j,ii) N2 N3 N4 N5];
                SSIM_RES = [SSIM_RES;err_threshold(ii), mod_threshold(jj) S1(i,j,ii) S2 S3 S4 S5];
                NMSE_RES = [NMSE_RES;err_threshold(ii), mod_threshold(jj) E1(i,j,ii) E2 E3 E4 E5];
                BPP_RES  = [BPP_RES;err_threshold(ii), mod_threshold(jj), ori_bpp, BPP1(i,j,ii), BPP2, BPP3, BPP4, BPP5];
                if alpha_1_3 + alpha_2_3 + alpha_1_4 + alpha_2_4 + alpha_1_5 + alpha_2_5 == 0
                else
                    Alphas = [Alphas; err_threshold(ii), mod_threshold(jj) alpha_1_3 alpha_2_3 alpha_1_4 alpha_2_4 alpha_1_5 alpha_2_5];
                end
                %ALL_STHRE_MTHRE_NSEBT = [ALL_STHRE_MTHRE_NSEBT; err_threshold(ii), mod_threshold(jj), N1(i,j,ii),S1(i,j,ii),E1(i,j,ii),ori_bpp,BPP1(i,j,ii),T1(i,j,ii),N2,S2,E2,ori_bpp,BPP2,T2];
                index = index + 1;
                k = k + 1;
            end
            
            
        end
        
        
        
        [OMP_ALG1, OMP_ALG2] = AAA_getOptDenoised_OneStage_PAP7_(im2uint8(X_OMP), im2uint8(Xhat_OMP));
        [PKbM_ALG1, PKbM_ALG2] = AAA_getOptDenoised_OneStage_PAP7_(im2uint8(X_OMP), im2uint8(X_PKbM));
        [PPKbM_1_ALG1, PPKbM_1_ALG2] = AAA_getOptDenoised_OneStage_PAP7_(im2uint8(X_OMP), im2uint8(X_PPKbM_1));
        [PPKbM_2_ALG1, PPKbM_2_ALG2] = AAA_getOptDenoised_OneStage_PAP7_(im2uint8(X_OMP), im2uint8(X_PPKbM_2));
        [UDbM_ALG1, UDbM_ALG2] = AAA_getOptDenoised_OneStage_PAP7_(im2uint8(X_OMP), im2uint8(X_UDbM));
        
        SSIM_OMP = [ssim_(Xhat_OMP,X_OMP), ssim_(OMP_ALG1,X_OMP),ssim_(OMP_ALG2,X_OMP)];
        SSIM_PKbM = [ssim_(X_PKbM,X_OMP), ssim_(PKbM_ALG1,X_OMP),ssim_(PKbM_ALG2,X_OMP)];
        SSIM_PPKbM_1 = [ssim_(X_PPKbM_1,X_OMP), ssim_(PPKbM_1_ALG1,X_OMP),ssim_(PPKbM_1_ALG2,X_OMP)];
        SSIM_PPKbM_2 = [ssim_(X_PPKbM_2,X_OMP), ssim_(PPKbM_2_ALG1,X_OMP),ssim_(PPKbM_2_ALG2,X_OMP)];
        SSIM_UDbM = [ssim_(X_UDbM,X_OMP), ssim_(UDbM_ALG1,X_OMP),ssim_(UDbM_ALG2,X_OMP)];
        
        PSNR_OMP = [psnr(im2uint8(Xhat_OMP),im2uint8(X_OMP)), psnr(im2uint8(OMP_ALG1),im2uint8(X_OMP)),psnr(im2uint8(OMP_ALG2),im2uint8(X_OMP))];
        PSNR_PKbM = [psnr(im2uint8(X_PKbM),im2uint8(X_OMP)), psnr(im2uint8(PKbM_ALG1),im2uint8(X_OMP)),psnr(im2uint8(PKbM_ALG2),im2uint8(X_OMP))];
        PSNR_PPKbM_1 = [psnr(im2uint8(X_PPKbM_1),im2uint8(X_OMP)), psnr(im2uint8(PPKbM_1_ALG1),im2uint8(X_OMP)),psnr(im2uint8(PPKbM_1_ALG2),im2uint8(X_OMP))];
        PSNR_PPKbM_2 = [psnr(im2uint8(X_PPKbM_2),im2uint8(X_OMP)), psnr(im2uint8(PPKbM_2_ALG1),im2uint8(X_OMP)),psnr(im2uint8(PPKbM_2_ALG2),im2uint8(X_OMP))];
        PSNR_UDbM = [psnr(im2uint8(X_UDbM),im2uint8(X_OMP)), psnr(im2uint8(UDbM_ALG1),im2uint8(X_OMP)),psnr(im2uint8(UDbM_ALG2),im2uint8(X_OMP))];
        
        NMSE_OMP = [immse(Xhat_OMP,X_OMP), immse(OMP_ALG1,X_OMP),immse(OMP_ALG2,X_OMP)];
        NMSE_PKbM = [immse(X_PKbM,X_OMP), immse(PKbM_ALG1,X_OMP),immse(PKbM_ALG2,X_OMP)];
        NMSE_PPKbM_1 = [immse(X_PPKbM_1,X_OMP), immse(PPKbM_1_ALG1,X_OMP),immse(PPKbM_1_ALG2,X_OMP)];
        NMSE_PPKbM_2 = [immse(X_PPKbM_2,X_OMP), immse(PPKbM_2_ALG1,X_OMP),immse(PPKbM_2_ALG2,X_OMP)];
        NMSE_UDbM = [immse(X_UDbM,X_OMP), immse(UDbM_ALG1,X_OMP),immse(UDbM_ALG2,X_OMP)];
        
        
        
        [~,Q_OMP_COEFF_ind] = imquantize(OMP_COEFF,-5:0.0196:5,0:1:511);
        [~,Q_PKbM_COEFF_ind] = imquantize(PKbM_COEFF,-5:0.0196:5,0:1:511);
        [~,Q_PPKbM_1_COEFF_ind] = imquantize(PPKbM_1_COEFF,-5:0.0196:5,0:1:511);
        [~,Q_PPKbM_2_COEFF_ind] = imquantize(PPKbM_2_COEFF,-5:0.0196:5,0:1:511);
        [~,Q_UDbM_COEFF_ind] = imquantize(UDbM_COEFF,-5:0.0196:5,0:1:511);
        
        [~, ~, Q_OMP_COEFF_bpp ] = get_prob_and_entropy(Q_OMP_COEFF_ind);
        [~, ~, Q_PKbM_COEFF_bpp ] = get_prob_and_entropy(Q_PKbM_COEFF_ind);
        [~, ~, Q_PPKbM_1_COEFF_bpp ] = get_prob_and_entropy(Q_PPKbM_1_COEFF_ind);
        [~, ~, Q_PPKbM_2_COEFF_bpp ] = get_prob_and_entropy(Q_PPKbM_2_COEFF_ind);
        [~, ~, Q_UDbM_COEFF_bpp ] = get_prob_and_entropy(Q_UDbM_COEFF_ind);
        
        [~, ~, OMP_IND_bpp ] = get_prob_and_entropy(OMP_IND);
        [~, ~, PKbM_IND_bpp ] = get_prob_and_entropy(PKbM_IND);
        [~, ~, PPKbM_1_IND_bpp ] = get_prob_and_entropy(PPKbM_1_IND);
        [~, ~, PPKbM_2_IND_bpp ] = get_prob_and_entropy(PPKbM_2_IND);
        [~, ~, UDbM_IND_bpp ] = get_prob_and_entropy(UDbM_IND);
        
        
    
        BPP_AVG = sum(BPP_RES(end-256+1:end,:),1)/256;
        BPP_ENT_AVG = [BPP_AVG(1,[1 2]), Q_OMP_COEFF_bpp+OMP_IND_bpp, Q_PKbM_COEFF_bpp+PKbM_IND_bpp, Q_PPKbM_1_COEFF_bpp+PPKbM_1_IND_bpp,Q_PPKbM_2_COEFF_bpp+PPKbM_2_IND_bpp,Q_UDbM_COEFF_bpp+UDbM_IND_bpp  ];
        TIME_AVG = sum(TIME_RES(end-256+1:end,:),1)/256;
        PARA_AVG = sum(PARA_RES(end-256+1:end,:),1)/256;
        SSIM_AVG = sum(SSIM_RES(end-256+1:end,:),1)/256;
        NMSE_AVG = sum(NMSE_RES(end-256+1:end,:),1)/256;
        ALPHAS_AVG = sum(Alphas(1:end,:),1)/size(Alphas,1);
        
        OVERALL_TIME = [OVERALL_TIME; TIME_AVG];
        OVERALL_PARA = [OVERALL_PARA; PARA_AVG];
        OVERALL_SSIM_avg = [OVERALL_SSIM_avg; SSIM_AVG];
        OVERALL_SSIM_FIN = [OVERALL_SSIM_FIN; SSIM_OMP SSIM_PKbM SSIM_PPKbM_1 SSIM_PPKbM_2 SSIM_UDbM];        
        OVERALL_NMSE_avg = [OVERALL_NMSE_avg; NMSE_AVG];
        OVERALL_NMSE_FIN = [OVERALL_NMSE_FIN; NMSE_OMP NMSE_PKbM NMSE_PPKbM_1 NMSE_PPKbM_2 NMSE_UDbM];
        OVERALL_PSNR_FIN = [OVERALL_PSNR_FIN; PSNR_OMP PSNR_PKbM PSNR_PPKbM_1 PSNR_PPKbM_2 PSNR_UDbM];
        OVERALL_BPP  = [OVERALL_BPP; BPP_AVG, BPP_ENT_AVG];
        OVERALL_ALPHAS = [OVERALL_ALPHAS; ALPHAS_AVG];
        Alphas = [];
        
        
        do_one_time = 1;
        
        
        
        check = 0;
    end
end
fprintf(strcat(img,' ... (Ours Alg. Reached %.2f%%)\n'),30*100/30);%2012

RESULTS.OVERALL_TIME = OVERALL_TIME;
RESULTS.OVERALL_PARA = OVERALL_PARA;
RESULTS.OVERALL_SSIM_avg = OVERALL_SSIM_avg;
RESULTS.OVERALL_SSIM_FIN = OVERALL_SSIM_FIN;
RESULTS.OVERALL_NMSE_avg = OVERALL_NMSE_avg;
RESULTS.OVERALL_NMSE_FIN = OVERALL_NMSE_FIN;
RESULTS.OVERALL_BPP  = OVERALL_BPP;
RESULTS.OVERALL_ALPHAS = OVERALL_ALPHAS;
RESULTS.JPEG_RES = JPEG_RES;
RESULTS.JPEG2000_RES = JPEG2000_RES;
RESULTS.ori_bpp = ori_bpp;
check = 0;

end


function       [entropy_, total_bits, bpp_entropy ] = get_prob_and_entropy(sequence)
sequence_ = sequence;
output_ = [];
while isempty(sequence_) == 0
    ind = find(sequence_(:) == sequence_(1));
    output_ = [output_;sequence_(1) size(ind,1) size(ind,1)/size(sequence,1)];
    sequence_(ind) = [];
    check = 0;
end
entropy_ = -1 * sum(output_(:,3).* log2(output_(:,3)));
total_bits = entropy_ * size(sequence,1);%size(output_,1);
bpp_entropy = total_bits /(256*256);
check = 0;
end
function [VecMatrix,kr,kc ] = AAA_getIMGvectors_(IMG,batchWidth)
VecLength = batchWidth;
RowMax = size(IMG,1)/VecLength;
ColMax = size(IMG,2);

kr = floor(RowMax);
kc = floor(ColMax);

for i=1:1:kr
    for j=1:1:kc
        VecMatrix(i,j,:)= IMG((i - 1) * VecLength + 1 : i * VecLength,j);
    end
    
end

end
function [RESULTS_fourth, coeff, iopt,  OMP_curERR, omp_ssim, y_i_omp] = AAA_SSIM_OMP_new_pap7_(A, y,k,ssim_threshold,err_threshold)

RESULTS_fourth = struct('YFIT',[],'R',[], ...
    'COEFF',[],'IOPT',[],'EOMP',[],'errSTRUCT',[],'TOMP',[],...
    'SNROMP',[],'curERR',[],'exact_sparsity',[],'qVals',[],...
    'qVals8',[],'qCodes',[],'newCOEFFs',[],'local_ssim',[],'local_ssim8',[]);

threshold_code = 1; %sparsity
switch threshold_code
    case 0
        stop = 0.002;
    case 1
        stop = k;
    case 2
        stop = k; %ssim_threshold;
end

xbeg = zeros(size(A,2),1);
eli = zeros(size(A,1),1);
B = A;
support=[];
temp=y;
normy = norm(y);
count = 1;
curr_level = 0;
err_ = 100;
local_ssim_old(count,1) = 0;
while curr_level < stop && err_>err_threshold %local_ssim_old(end,1)<= ssim_threshold
    % tic;
    
    ST = abs(B' * temp);
    [a, b] = max(ST);
    support = [support b];
    xfinal = A(:, support)\y;
    
    
    y_hat = A(:,support) * xfinal;
    if length(xfinal)==1
        y_hat_temp = y_hat;
    end
    temp = y - y_hat;
    %%%%%%%%%
    curERR(count,1) = norm(temp);
    err_ = curERR(count,1)*100/normy;
    %local_ssim_old(count,1) = ssim_(im2uint8(y_hat),im2uint8(y));
    
    y_i_omp(:,count) = y_hat;
    if length(xfinal)>=1
        
        OMP_curERR(count,1) = immse(y_hat,y);
        
    else
        
    end
    
    
    switch threshold_code
        case 0
            curr_level = curERR(count,1);
        case 1
            curr_level = count;
        case 2
            curr_level = local_ssim(count,1);
    end
    count = count + 1;
    
    
end

coeff = xfinal;

iopt = support';
omp_ssim = ssim_(y_hat,y);%local_ssim_old;

RESULTS_fourth.YFIT = A(:,support) * xfinal;
RESULTS_fourth.y_i_omp = y_i_omp;

RESULTS_fourth.R = y - RESULTS_fourth.YFIT;
RESULTS_fourth.SNROMP = snr(y,RESULTS_fourth.R);
RESULTS_fourth.IOPT = support';
RESULTS_fourth.COEFF = xfinal;
RESULTS_fourth.exact_sparsity = length(support);
RESULTS_fourth.curERR = curERR;

RESULTS_fourth.OMP_curERR= OMP_curERR;


temp1 = omp_ssim;
omp_ssim = [];
omp_ssim = temp1;

end
function [RESULTS_fourth, coeff, iopt,  OMP_curERR, omp_ssim, y_i_omp] = Enh_Convergence_OMP_pap7_(A, y,k,y_init,ssim_threshold,err_threshold)

RESULTS_fourth = struct('YFIT',[],'R',[], ...
    'COEFF',[],'IOPT',[],'EOMP',[],'errSTRUCT',[],'TOMP',[],...
    'SNROMP',[],'curERR',[],'exact_sparsity',[],'qVals',[],...
    'qVals8',[],'qCodes',[],'newCOEFFs',[],'local_ssim',[],'local_ssim8',[]);


threshold_code = 1; %sparsity
switch threshold_code
    case 0
        stop = 0.002;
    case 1
        stop = k;
    case 2
        stop = k; %ssim_threshold;
end

xbeg = zeros(size(A,2),1);
eli = zeros(size(A,1),1);
B = A;
support=[];
temp=y-y_init;
normy = norm(y);
err_ = 100;
count = 1;
curr_level = 0;
[AA, BB] = size(A);
R = zeros(AA,BB);
local_ssim_old(count,1) = 0;
while curr_level < stop && err_>err_threshold
    ST = abs(B' * temp);
    [a, b] = max(ST);
    support = [support b];
    xfinal = A(:, support)\(y-y_init);
    e_hat = A(:,support) * xfinal;
    e_hat2 = A(:,support) * xfinal;
    if length(xfinal)==1
        y_hat_temp = e_hat;
    end
    temp = y - (e_hat +y_init) ;
    curERR(count,1) = norm(temp);
    err_ = curERR(count,1)*100/normy;
    %local_ssim_old(count,1) = ssim_(im2uint8(e_hat + y_init),im2uint8(y));
    
    y_i_omp(:,count) = e_hat + y_init;
    
    if length(xfinal)>=1
        
        OMP_curERR(count,1) = immse(e_hat+y_init,y);
        
    else
        
    end
    
    switch threshold_code
        case 0
            curr_level = curERR(count,1);
        case 1
            curr_level = count;
        case 2
            curr_level = local_ssim(count,1);
    end
    count = count + 1;
    
end

coeff = xfinal;

iopt = support';
omp_ssim = ssim_(e_hat + y_init,y);
RESULTS_fourth.YFIT = A(:,support) * xfinal;
RESULTS_fourth.y_i_omp = y_i_omp;
RESULTS_fourth.R = y - RESULTS_fourth.YFIT;
RESULTS_fourth.SNROMP = snr(y,RESULTS_fourth.R);
RESULTS_fourth.IOPT = support';
RESULTS_fourth.COEFF = xfinal;
RESULTS_fourth.exact_sparsity = length(support);
RESULTS_fourth.curERR = curERR;

temp1 = omp_ssim;
omp_ssim = [];
omp_ssim = temp1;

end
function [BM3D_, TVL1_] = AAA_getOptDenoised_OneStage_PAP7_(origIm, noisedIm)
SSIM_ = [];
NMSE_ = [];
Error = noisedIm - origIm;
stdev = std2(Error);
tic;
[v1, BM3D_] = BM3D(1,noisedIm,stdev,'np',0);

TVL1_ = TVL1denoise(noisedIm, 2,20);
check = 0;

end
function [Dic] = getDic_(dicType,vecLength)
if dicType == 'LDic'
    
    dd = load('MOD_LEARNED_256_512.mat');
    Dic = dd.Dic;
else
    
    [Dic_r,nbvect] = wmpdictionary(vecLength,'lstcpt',{{'wphaar',2},'dct',});%{'wphaar',2},,{'wpsym4',2}}'dct'
    aa = ones(vecLength,1);
    Dic = full(Dic_r); %[aa_n full(Dic_r)];
    
end
end
function Dic = getLearnedAtoms_()
img1 = imread('(250012).jpg');%imread('optImage_trainSET.tif');
img_double = im2double(img1);
imgDouble(:,:)=img_double(:,:,1);
[Batches, batchPerRow,batchPerCol] = getBatches_(imgDouble,16);
k = 1;
for i=1:1:batchPerRow
    for j=1:1:batchPerCol
        Vector1(:,k) = batchToVector_(Batches(i,j,:,:),16);
        k = k + 1;
    end
end
img2 = imread('(250014).jpg');%imread('optImage_trainSET.tif');
img_double = im2double(img2);
imgDouble(:,:)=img_double(:,:,1);
[Batches, batchPerRow,batchPerCol] = getBatches_(imgDouble,16);
k = 1;
for i=1:1:batchPerRow
    for j=1:1:batchPerCol
        Vector2(:,k) = batchToVector_(Batches(i,j,:,:),16);
        k = k + 1;
    end
end

img3 = imread('(250005).jpg');%imread('optImage_trainSET.tif');
img_double = im2double(img3);
imgDouble(:,:)=img_double(:,:,1);
[Batches, batchPerRow,batchPerCol] = getBatches_(imgDouble,16);
k = 1;
for i=1:1:batchPerRow
    for j=1:1:batchPerCol
        Vector3(:,k) = batchToVector_(Batches(i,j,:,:),16);
        k = k + 1;
    end
end

img4 = imread('(250036).jpg');%imread('optImage_trainSET.tif');
img_double = im2double(img4);
imgDouble(:,:)=img_double(:,:,1);
[Batches, batchPerRow,batchPerCol] = getBatches_(imgDouble,16);
k = 1;
for i=1:1:batchPerRow
    for j=1:1:batchPerCol
        Vector4(:,k) = batchToVector_(Batches(i,j,:,:),16);
        k = k + 1;
    end
end
img5 = imread('(250042).jpg');%imread('optImage_trainSET.tif');
img_double = im2double(img5);
imgDouble(:,:)=img_double(:,:,1);
[Batches, batchPerRow,batchPerCol] = getBatches_(imgDouble,16);
k = 1;
for i=1:1:batchPerRow
    for j=1:1:batchPerCol
        Vector5(:,k) = batchToVector_(Batches(i,j,:,:),16);
        k = k + 1;
    end
end

k = 1;
for i=1:1:batchPerRow
    for j=1:1:batchPerCol
        Vector6(:,k) = batchToVector_(Batches(i,j,:,:),16);
        k = k + 1;
    end
end


Vector = [Vector1 Vector2 Vector3 Vector4 Vector5];


options.niter_learning = 20;
options.K = 512;
options.nbr_max_atoms = 12;
[Dic,X,E] = perform_dictionary_learning(Vector, options);
x=0;
save('MOD_LEARNED_256','Dic');
y = 0;


end
function [Vector] = batchToVector_(Batch, batchSize)%col to vectors
vector = [];
for row=1:1:batchSize
    for col = 1:1:batchSize
        vector = [vector Batch(1,1,row,col)];
    end
    
end
Vector = vector';

end
function [batchMatrix, batchPerRow,batchPerCol] = getBatches_(singleColorLevelMatrix,batchSize)

batchPerRow = size(singleColorLevelMatrix,1)/batchSize;
batchPerCol = size(singleColorLevelMatrix,2)/batchSize;

kr = 0:1:batchPerRow;
kc = 0:1:batchPerCol;
for i=1:1:floor(batchPerRow)
    for j=1:1:floor(batchPerCol)
        batchMatrix(i,j,:,:)= singleColorLevelMatrix(kr(i)*batchSize+1:kr(i)*batchSize+batchSize,kc(j)*batchSize+1:kc(j)*batchSize+batchSize);
    end
    
end

end
function [ssimval, ssimmap] = ssim_(varargin)
%SSIM Structural Similarity Index for measuring image quality
%   SSIMVAL = SSIM(A, REF) calculates the Structural Similarity Index
%   (SSIM) value for image A, with the image REF as the reference. A and
%   REF can be 2D grayscale or 3D volume images, and must be of the same
%   size and class.
%
%   [SSIMVAL, SSIMMAP] = SSIM(A, REF) also returns the local SSIM value for
%   each pixel in SSIMMAP. SSIMMAP has the same size as A.
%
%   [SSIMVAL, SSIMMAP] = SSIM(A, REF, NAME1, VAL1,...) calculates the SSIM
%   value using name-value pairs to control aspects of the computation.
%   Parameter names can be abbreviated.
%
%   Parameters include:
%
%   'Radius'                 - Specifies the standard deviation of
%                              isotropic Gaussian function used for
%                              weighting the neighborhood pixels around a
%                              pixel for estimating local statistics. This
%                              weighting is used to avoid blocking
%                              artifacts in estimating local statistics.
%                              The default value is 1.5.
%
%   'DynamicRange'           - Positive scalar, L, that specifies the
%                              dynamic range of the input image. By
%                              default, L is chosen based on the class of
%                              the input image A, as L =
%                              diff(getrangefromclass(A)). Note that when
%                              class of A is single or double, L = 1 by
%                              default.
%
%   'RegularizationConstants'- Three-element vector, [C1 C2 C3], of
%                              non-negative real numbers that specifies the
%                              regularization constants for the luminance,
%                              contrast, and structural terms (see [1]),
%                              respectively. The regularization constants
%                              are used to avoid instability for image
%                              regions where the local mean or standard
%                              deviation is close to zero. Therefore, small
%                              non-zero values should be used for these
%                              constants. By default, C1 = (0.01*L).^2, C2
%                              = (0.03*L).^2, and C3 = C2/2, where L is the
%                              specified 'DynamicRange' value. If a value
%                              of 'DynamicRange' is not specified, the
%                              default value is used (see name-value pair
%                              'DynamicRange').
%
%   'Exponents'               - Three-element vector [alpha beta gamma],
%                               of non-negative real numbers that specifies
%                               the exponents for the luminance, contrast,
%                               and structural terms (see [1]),
%                               respectively. By default, all the three
%                               exponents are 1, i.e. the vector is [1 1
%                               1].
%
%   Notes
%   -----
%   1. A and REF can be arrays of upto three dimensions. All 3D arrays
%      are considered 3D volumetric images. RGB images will also be
%      processed as 3D volumetric images.
%
%   2. Input image A and reference image REF are converted to
%      floating-point type for internal computation.
%
%   3. For signed-integer images (int16), an offset is applied to bring the
%      gray values in the non-negative range before computing the SSIM
%      index.
%
%   Example
%   ---------
%   This example shows how to compute SSIM value for a blurred image given
%   the original reference image.
%
%   ref = imread('pout.tif');
%   A = imgaussfilt(ref, 1.5, 'FilterSize', 11, 'Padding', 'replicate');
%
%   subplot(1,2,1); imshow(ref); title('Reference Image');
%   subplot(1,2,2); imshow(A);   title('Blurred Image');
%
%   [ssimval, ssimmap] = ssim(A,ref);
%
%   fprintf('The SSIM value is %0.4f.\n',ssimval);
%
%   figure, imshow(ssimmap,[]);
%   title(sprintf('SSIM Index Map - Mean SSIM Value is %0.4f',ssimval));
%
%   Class Support
%   -------------
%   Input arrays A and REF must be one of the following classes: uint8,
%   int16, uint16, single, or double. Both A and REF must be of the same
%   class. They must be nonsparse. SSIMVAL is a scalar and SSIMMAP is an
%   array of the same size as A. Both SSIMVAL and SSIMMAP are of class
%   double, unless A and REF are of class single in which case SSIMVAL and
%   SSIMMAP are of class single.
%
%   References:
%   -----------
%   [1] Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
%       Quality Assessment: From Error Visibility to Structural
%       Similarity," IEEE Transactions on Image Processing, Volume 13,
%       Issue 4, pp. 600- 612, 2004.
%
%   See also IMMSE, MEAN, MEDIAN, PSNR, SUM, VAR.

%   Copyright 2013-2017 The MathWorks, Inc.

narginchk(2,10);

args = matlab.images.internal.stringToChar(varargin);
[A, ref, C, exponents, radius] = parse_inputs(args{:});

if isempty(A)
    ssimval = zeros(0, 'like', A);
    ssimmap = A;
    return;
end

if isa(A,'int16') % int16 is the only allowed signed-integer type for A and ref.
    % Add offset for signed-integer types to bring values in the
    % non-negative range.
    A = double(A) - double(intmin('int16'));
    ref = double(ref) - double(intmin('int16'));
elseif isinteger(A)
    A = double(A);
    ref = double(ref);
end

filtRadius = ceil(radius*3); % 3 Standard deviations include >99% of the area.
filtSize = 2*filtRadius + 1;

if ismatrix(A)
    gaussFilterFcn = @(X)imgaussfilt(X, radius, 'FilterSize', filtSize, 'Padding','replicate');
else
    gaussFilterFcn = @(X)imgaussfilt3(X, radius, 'FilterSize', filtSize, 'Padding','replicate');
end
% Weighted-mean and weighted-variance computations
mux2 = gaussFilterFcn(A);
muy2 = gaussFilterFcn(ref);
muxy = mux2.*muy2;
mux2 = mux2.^2;
muy2 = muy2.^2;

sigmax2 = gaussFilterFcn(A.^2) - mux2;
sigmay2 = gaussFilterFcn(ref.^2) - muy2;
sigmaxy = gaussFilterFcn(A.*ref) - muxy;

% Compute SSIM index
if (C(3) == C(2)/2) && isequal(exponents(:),ones(3,1))
    % Special case: Equation 13 from [1]
    num = (2*muxy + C(1)).*(2*sigmaxy + C(2));
    den = (mux2 + muy2 + C(1)).*(sigmax2 + sigmay2 + C(2));
    if (C(1) > 0) && (C(2) > 0)
        ssimmap = num./den;
    else
        % Need to guard against divide-by-zero if either C(1) or C(2) is 0.
        isDenNonZero = (den ~= 0);
        ssimmap = ones(size(A));
        ssimmap(isDenNonZero) = num(isDenNonZero)./den(isDenNonZero);
    end
    
else
    % General case: Equation 12 from [1]
    % Luminance term
    if (exponents(1) > 0)
        num = 2*muxy + C(1);
        den = mux2 + muy2 + C(1);
        ssimmap = guardedDivideAndExponent(num,den,C(1),exponents(1));
    else
        ssimmap = ones(size(A), 'like', A);
    end
    
    % Contrast term
    sigmaxsigmay = [];
    if (exponents(2) > 0)
        sigmaxsigmay = sqrt(sigmax2.*sigmay2);
        num = 2*sigmaxsigmay + C(2);
        den = sigmax2 + sigmay2 + C(2);
        ssimmap = ssimmap.*guardedDivideAndExponent(num,den,C(2),exponents(2));
    end
    
    % Structure term
    if (exponents(3) > 0)
        num = sigmaxy + C(3);
        if isempty(sigmaxsigmay)
            sigmaxsigmay = sqrt(sigmax2.*sigmay2);
        end
        den = sigmaxsigmay + C(3);
        ssimmap = ssimmap.*guardedDivideAndExponent(num,den,C(3),exponents(3));
    end
    
end

ssimval = mean(ssimmap(:));

end
function component = guardedDivideAndExponent(num, den, C, exponent)

if C > 0
    component = num./den;
else
    component = ones(size(num),'like',num);
    isDenNonZero = (den ~= 0);
    component(isDenNonZero) = num(isDenNonZero)./den(isDenNonZero);
end

if (exponent ~= 1)
    component = component.^exponent;
end

end
function [A, ref, C, exponents, radius] = parse_inputs(varargin)

validImageTypes = {'uint8','uint16','int16','single','double'};

A = varargin{1};
validateattributes(A,validImageTypes,{'nonsparse','real'},mfilename,'A',1);

ref = varargin{2};
validateattributes(ref,validImageTypes,{'nonsparse','real'},mfilename,'REF',2);

if ~isa(A,class(ref))
    error(message('images:validate:differentClassMatrices','A','REF'));
end

if ~isequal(size(A),size(ref))
    error(message('images:validate:unequalSizeMatrices','A','REF'));
end

if (ndims(A) > 3)
    error(message('images:validate:tooManyDimensions','A and REF',3));
end

% Default values for parameters
dynmRange = diff(getrangefromclass(A));
C = [];
exponents = [1 1 1];
radius = 1.5;

args_names = {'dynamicrange', 'regularizationconstants','exponents',...
    'radius'};

for i = 3:2:nargin
    arg = varargin{i};
    if ischar(arg)
        idx = find(strncmpi(arg, args_names, numel(arg)));
        if isempty(idx)
            error(message('images:validate:unknownInputString', arg))
            
        elseif numel(idx) > 1
            error(message('images:validate:ambiguousInputString', arg))
            
        elseif numel(idx) == 1
            if (i+1 > nargin)
                error(message('images:validate:missingParameterValue'));
            end
            if idx == 1
                dynmRange = varargin{i+1};
                validateattributes(dynmRange,{'numeric'},{'positive', ...
                    'finite', 'real', 'nonempty','scalar'}, mfilename, ...
                    'DynamicRange',i);
                dynmRange = double(dynmRange);
                
            elseif idx == 2
                C = varargin{i+1};
                validateattributes(C,{'numeric'},{'nonnegative','finite', ...
                    'real','nonempty','vector', 'numel', 3}, mfilename, ...
                    'RegularizationConstants',i);
                C = double(C);
                
            elseif idx == 3
                exponents = varargin{i+1};
                validateattributes(exponents,{'numeric'},{'nonnegative', ...
                    'finite', 'real', 'nonempty','vector', 'numel', 3}, ...
                    mfilename,'Exponents',i);
                exponents = double(exponents);
                
            elseif idx == 4
                radius = varargin{i+1};
                validateattributes(radius,{'numeric'},{'positive','finite', ...
                    'real', 'nonempty','scalar'}, mfilename,'Radius',i);
                radius = double(radius);
            end
        end
    else
        error(message('images:validate:mustBeString'));
    end
end

% If 'RegularizationConstants' is not specified, choose default C.
if isempty(C)
    C = [(0.01*dynmRange).^2 (0.03*dynmRange).^2 ((0.03*dynmRange).^2)/2];
end

end
function [ssimval, ssimmap] = ssim(varargin)
%SSIM Structural Similarity Index for measuring image quality
%   SSIMVAL = SSIM(A, REF) calculates the Structural Similarity Index
%   (SSIM) value for image A, with the image REF as the reference. A and
%   REF can be 2D grayscale or 3D volume images, and must be of the same
%   size and class. The similarity metric, SSIMVAL, is a double valued
%   scalar. A value of 1 corresponds to the highest quality (when A and
%   REF are equivalent).
%
%   [SSIMVAL, SSIMMAP] = SSIM(A, REF) also returns the local SSIM value for
%   each pixel in SSIMMAP. SSIMMAP has the same size as A.
%
%   [SSIMVAL, SSIMMAP] = SSIM(A, REF, NAME1, VAL1,...) calculates the SSIM
%   value using name-value pairs to control aspects of the computation.
%   Parameter names can be abbreviated.
%
%   Parameters include:
%
%   'Radius'                 - Specifies the standard deviation of
%                              isotropic Gaussian function used for
%                              weighting the neighborhood pixels around a
%                              pixel for estimating local statistics. This
%                              weighting is used to avoid blocking
%                              artifacts in estimating local statistics.
%                              The default value is 1.5.
%
%   'DynamicRange'           - Positive scalar, L, that specifies the
%                              dynamic range of the input image. By
%                              default, L is chosen based on the class of
%                              the input image A, as L =
%                              diff(getrangefromclass(A)). Note that when
%                              class of A is single or double, L = 1 by
%                              default.
%
%   'RegularizationConstants'- Three-element vector, [C1 C2 C3], of
%                              non-negative real numbers that specifies the
%                              regularization constants for the luminance,
%                              contrast, and structural terms (see [1]),
%                              respectively. The regularization constants
%                              are used to avoid instability for image
%                              regions where the local mean or standard
%                              deviation is close to zero. Therefore, small
%                              non-zero values should be used for these
%                              constants. By default, C1 = (0.01*L).^2, C2
%                              = (0.03*L).^2, and C3 = C2/2, where L is the
%                              specified 'DynamicRange' value. If a value
%                              of 'DynamicRange' is not specified, the
%                              default value is used (see name-value pair
%                              'DynamicRange').
%
%   'Exponents'               - Three-element vector [alpha beta gamma],
%                               of non-negative real numbers that specifies
%                               the exponents for the luminance, contrast,
%                               and structural terms (see [1]),
%                               respectively. By default, all the three
%                               exponents are 1, i.e. the vector is [1 1
%                               1].
%
%   'DataFormat'              - Dimension labels of the input data A and REF
%                               specified as a string scalar or character
%                               vector. The format options 'S','C', and 'B'
%                               are supported. The options 'S', 'C' and 'B'
%                               correspond to spatial, channel, and batch
%                               dimensions, respectively. A separate SSIMVA
%                               SSIMVAL and SSIMMAP output will be returned
%                               for each non-spatial dimension.
%
%   Class Support
%   -------------
%   Input arrays A and REF must be one of the following classes: uint8,
%   int16, uint16, single, or double. Both A and REF must be of the same
%   class. They must be nonsparse. SSIMVAL is a scalar and SSIMMAP is an
%   array of the same size as A. Both SSIMVAL and SSIMMAP are of class
%   double, unless A and REF are of class single in which case SSIMVAL and
%   SSIMMAP are of class single.
%
%   Notes
%   -----
%   1. A and REF can be arrays of up to three dimensions. All 3D arrays
%      are considered 3D volumetric images. RGB images will also be
%      processed as 3D volumetric images.
%
%   2. Input image A and reference image REF are converted to
%      floating-point type for internal computation.
%
%   3. For signed-integer images (int16), an offset is applied to bring the
%      gray values in the non-negative range before computing the SSIM
%      index.
%
%   4. SSIMVAL is 1 when A and REF are equivalent indicating highest
%      quality. Smaller values indicate deviations. For some combinations
%      of inputs and parameters, SSIMVAL can be negative.
%
%   5. When non-integer valued 'Exponents' are used, intermediate
%      luminance, contrast and structural terms are clamped to [0, inf]
%      range to prevent complex valued outputs.
%
%   References:
%   -----------
%   [1] Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
%       Quality Assessment: From Error Visibility to Structural
%       Similarity," IEEE Transactions on Image Processing, Volume 13,
%       Issue 4, pp. 600- 612, 2004.
%
%   Example
%   ---------
%   % This example shows how to compute SSIM value for a blurred image
%   % given the original reference image.
%
%   ref = imread('pout.tif');
%   A = imgaussfilt(ref, 1.5, 'FilterSize', 11, 'Padding', 'replicate');
%
%   subplot(1,2,1); imshow(ref); title('Reference Image');
%   subplot(1,2,2); imshow(A);   title('Blurred Image');
%
%   [ssimval, ssimmap] = ssim(A,ref);
%
%   fprintf('The SSIM value is %0.4f.\n',ssimval);
%
%   figure, imshow(ssimmap,[]);
%   title(sprintf('SSIM Index Map - Mean SSIM Value is %0.4f',ssimval));
%
%   See also IMMSE, MEAN, MEDIAN, PSNR, SUM, VAR.

%   Copyright 2013-2021 The MathWorks, Inc.

narginchk(2,12);

args = matlab.images.internal.stringToChar(varargin);

[A, ref, C, exponents, radius, filtSize, dataformat] = images.internal.qualitymetric.ssimParseInputs(args{:});

% Validate ndims of A differently depending on whether or not DataFormat
% specified to preserve behavior pre introduction of DataFormat Name/Value.
if isempty(dataformat)
    if (ndims(A) > 3)
        error(message('images:ssim:supportedDimsUnformattedInput'))
    else
        dataformat = repmat('S',1,ndims(A));
    end
else
    if (sum(dataformat == 'S') > 3) || (sum(dataformat == 'S') < 2)
        error(message('images:ssim:unsupportedDataFormatSpatialDims'));
    end
end

if isempty(A)
    ssimval = zeros(0, 'like', A);
    ssimmap = A;
    return;
end

numSpatialDims = sum(dataformat == 'S');
if numSpatialDims == 2
    gaussFilterFcn = @(X)imgaussfilt(X, radius, 'FilterSize', filtSize, 'Padding','replicate');
elseif numSpatialDims == 3
    gaussFilterFcn = @(X) spatialGaussianFilter(X, double(radius), double(filtSize), 'replicate');
else
    assert(false,'Unexpected number of spatial dimensions.');
end

[ssimval,ssimmap] = images.internal.qualitymetric.ssimalgo(A,ref,gaussFilterFcn,exponents,C,numSpatialDims);

end
function A = spatialGaussianFilter(A, sigma, hsize, padding)
[hCol,hRow,hSlc] = createSeparableGaussianKernel(sigma, hsize);

A = imfilter(A, hRow, padding, 'conv', 'same');
A = imfilter(A, hCol, padding, 'conv', 'same');
A = imfilter(A, hSlc, padding, 'conv', 'same');
end
function [hcol,hrow,hslc] = createSeparableGaussianKernel(sigma, hsize)

hcol = images.internal.createGaussianKernel(sigma(1), hsize(1));
hrow = reshape(hcol,1,[]);
hslc = reshape(hcol,1,1,[]);
end
function [alpha_vector, init] = prep_init_matrix(j,J_INDEX,X,X_HAT)

%fix_ind = [1:mod_:size(y_hat,2),size(y_hat,2)];
%init_mat = zeros(size(y_hat,1),size(y_hat,2));
%alpha_vector = zeros(size(y_hat,2),1);
%init_mat(:,fix_ind) = hat_mat(:,fix_ind);
for jj = 1:1:length(J_INDEX)-1
    if j>J_INDEX(jj) && j<J_INDEX(jj+1)
        i1 = J_INDEX(jj);
        i2 = J_INDEX(jj+1);
        y_i1 = X_HAT(:,j-1);%/norm(X_HAT(:,j-1),2);
        y_i2 = X_HAT(:,i2);%/norm(X_HAT(:,i2),2);
        alpha_vector = inv([y_i1, y_i2]'*[y_i1, y_i2])*[y_i1, y_i2]'* X ;%init_mat(:,j) = 1 * ((j-i1)/(i2-i1-1) * y_hat(:,i2) + (1-(j-i1)/(i2-i1-1)) * y_hat(:,i1));
        
        init =[y_i1, y_i2]*alpha_vector;
        check = 0;
        break;
    end
end



check = 0;
end
function [Atom_1,Atom_2] = get_2_bases(j,J_INDEX,X,X_HAT)
%fix_ind = [1:mod_:size(y_hat,2),size(y_hat,2)];
%init_mat = zeros(size(y_hat,1),size(y_hat,2));
%alpha_vector = zeros(size(y_hat,2),1);
%init_mat(:,fix_ind) = hat_mat(:,fix_ind);
%for j = 1:1:size(y_hat,2)
for jj = 1:1:length(J_INDEX)-1
    if j>J_INDEX(jj) && j<J_INDEX(jj+1)
        i1 = J_INDEX(jj);
        i2 = J_INDEX(jj+1);
        Atom_1 = X_HAT(:,j-1);
        Atom_2 = X_HAT(:,i2);
        break;
        
    end
end



%end
end
function AA = dcblker(A)
B = ones(length(A),1);
DC = B./norm(B,2);
AA = A - (DC'*A)*DC;
check = 0;
end


function [alpha_1, alpha_2, init] = getAlpha(j,J_INDEX,X,X_HAT)
%fix_ind = [1:mod_:size(y_hat,2),size(y_hat,2)];
%init_mat = zeros(size(y_hat,1),size(y_hat,2));
%alpha_vector = zeros(size(y_hat,2),1);
%init_mat(:,fix_ind) = hat_mat(:,fix_ind);

for jj = 1:1:length(J_INDEX)-1
    if j>J_INDEX(jj) && j<J_INDEX(jj+1)
        i1 = J_INDEX(jj);
        i2 = J_INDEX(jj+1);
        err = [];
        alpha = 0:0.01:1;
        for m = 1:1:length(alpha)
            error = (alpha(m)*X_HAT(:,j-1)+(1-alpha(m))*X_HAT(:,i2))-X ;%init_mat(:,j) = 1 * ((j-i1)/(i2-i1-1) * y_hat(:,i2) + (1-(j-i1)/(i2-i1-1)) * y_hat(:,i1));
            err = [err;norm(error,2)/norm(X,2)];
        end
        [a b] = min(err);
        init = alpha(b)*X_HAT(:,j-1)+(1-alpha(b))*X_HAT(:,i2);
        alpha_1 = alpha(b);
        alpha_2 = 1 - alpha_1;
        check = 0;
        break;
    end
end




check = 0;
end