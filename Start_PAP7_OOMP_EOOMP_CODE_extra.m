%imDim,imName,imExt,dicType,ssimDIF, enh
function Start_PAP7_OOMP_EOOMP_CODE_extra(file_index)
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
    %the_size=AAA_getSize(fName);
    pre_res  = strcat('PAP7_RESULTS_DS3/PAP7_OOMP_EOOMP');
    dicType_ = 'SDic';
    %[STRUCTURED_RESULTS_SA] = AAA_Second_BE_Final_SparseALGs_enh_PAP5(fName,batchSize,sparse,dicType_);
    [SDIC_RESULTS] = PAP7_main_code(fName,ori_bpp,batchSize,sparse,dicType_,JPEG_RES,JPEG2000_RES);
    dicType_ = 'LDic';
    fName = strcat('(',int2str(i),')',imExt);
    %[LEARNED_RESULTS_SA] = AAA_Second_BE_Final_SparseALGs_enh_PAP5(fName,batchSize,sparse,dicType_);
    [LDIC_RESULTS] = PAP7_main_code(fName,ori_bpp,batchSize,sparse,dicType_,JPEG_RES,JPEG2000_RES);
    
    
    save(strcat(pre_res,'_file_','(',int2str(i),')'),'SDIC_RESULTS','LDIC_RESULTS');
    check = 0;
    
end
catch 
    fid = fopen( strcat('Avoid_(',int2str(file_index),')_.txt'), 'wt' );

end
end


function [RESULTS] = PAP7_main_code(img,ori_bpp,batchSize,sparse_,dicType,JPEG_RES,JPEG2000_RES) %#ok<*INUSD>
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
iter = 0;
vecLength = batchSize;%batchSize*batchSize;
img_double = im2double(img_pixels);
Dic = getDic_(dicType,vecLength);
[VecMatrix,kr,kc] = AAA_getIMGvectors_(img_double,batchSize);
index = 1;
warning('off','all');

ssim_threshold = 0.4:0.1:0.9;
err_threshold = 5:5:50;
mod_threshold = 2:2:30;
OVERALL_RESULTS = [];
OVERALL_SSIM_RESULTS = [];
OVERALL_PSNR_RESULTS = [];
OVERALL_NMSE_RESULTS = [];
OVERALL_BPP_RESULTS = [];

ALL_STHRE_MTHRE_NSEBT = [];
indicator_ = 1;
DicCoh = Dic'*Dic;
%omp_ci_matrix = zeros(size(Dic,2),kr*kc,length(err_threshold));
%eomp_ci_matrix = zeros(size(Dic,2),kr*kc,length(err_threshold),length(mod_threshold));
for ii = 1:1:length(err_threshold)
    do_one_time = 0;
    kk = 1;
    for jj = 1:1:length(mod_threshold)
        if do_one_time == 0
            X_OMP = [];
            OMP_COEFF = [];
            OMP_IND = [];
            Xhat_OMP = [];
            X_OMP_ = [];
            Xhat_OMP_ = [];
            Xhat_OMP_BRe = [];
            Xhat_OMP_BRe_ = [];
            N1 = zeros(kr,kc,length(err_threshold));
            S1  = zeros(kr,kc,length(err_threshold));
            E1  = zeros(kr,kc,length(err_threshold));
            BPP1 = zeros(kr,kc,length(err_threshold));
            T1 = zeros(kr,kc,length(err_threshold));
        end
        EOMP_COEFF = [];
        EOMP_IND = [];
        X_eOMP = [];
        Xhat_eOMP = [];
        Xhat_eOMP_BRe  = [];
        X_eOMP_ = [];
        Xhat_eOMP_ = [];
        Xhat_eOMP_BRe_ = [];
    
        y_init = zeros(batchSize,1);
        clc;
        fprintf(strcat(img,' ... ',dicType, ' ... ', '(Ours Alg. Reached %.2f%%)\n'),(indicator_-1)*100/(length(err_threshold)*length(mod_threshold)));%2012
        indicator_ = indicator_ + 1;
        
        kkk = 1;
        for i=1:1:kr
            k = 1;
            for j=1:1:kc
                Vector =  im2double(img_pixels(256*(i-1)+1:i*256,j)); %VecMatrix(i,j,:);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(OMP(1990) - StOMP(2006)? - SWOMP(2009)  LASSO(2009) - COSAMP(2008) - GSM(2020)-
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OGA() - OMPvar(2015) - SAEN(2018) - GOMP(2010) -  - cpwwen(2018) - clarswlasso(2018) )%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if do_one_time == 0
                    tic;
                    %[~, omp_coeff, omp_iopt, omp_error,omp_ssim, omp_hat]       = AAA_SSIM_OMP_new_pap7_(Dic, Vector,sparse,1, err_threshold(ii));
                    %[~, ~, ~, ~,omp_ssim, ~]       = AAA_SSIM_OMP_new_pap7_(Dic, Vector,sparse,1, err_threshold(ii));
                    [~, ~, ~,~, omp_coeff,omp_iopt,~, omp_hat,omp_error,omp_ssim] = OOMPex_PAP7(Vector, Dic, err_threshold(ii), sparse);

                    if size(omp_iopt,1) == 1
                        SIZE_2 = (ceil(log2(size(omp_coeff,1))) + bits_per_coeff + (size(omp_iopt,1)-1)*(bits_per_coeff + bits_per_index));
                        Approx_ = omp_hat(:,end);
                    else
                        y_init_ = zeros(size(Vector,1),1);
                        [~,SIZE_2,Approx_,~] = runBRe(Vector,omp_error,omp_coeff,omp_iopt,Dic, Vector-omp_hat(:,end), size(omp_iopt,1), DicCoh,bits_per_index,bits_per_coeff,0,omp_hat(:,end-1),y_init_);
                    end
                    OMP_COEFF = [OMP_COEFF;omp_coeff];
                    OMP_IND = [OMP_IND;omp_iopt];
                    T1(i,j,ii) = toc;
                    %omp_ci_matrix(omp_iopt,kk,ii)=omp_coeff;
                    N1(i,j,ii) = size(omp_iopt,1);
                    S1(i,j,ii) = omp_ssim(end,1);
                    E1(i,j,ii) = omp_error(end,1);
                    BPP1(i,j,ii) = (ceil(log2(size(omp_coeff,1))) + bits_per_coeff + (N1(i,j,ii)-1)*(bits_per_coeff + bits_per_index))/256;
                    BPP1_BRe(i,j,ii) = SIZE_2/256;
                    X_OMP = [X_OMP Vector];
                    Xhat_OMP = [Xhat_OMP omp_hat(:,end)];
                    Xhat_OMP_BRe = [Xhat_OMP_BRe Approx_];
                    
                    kk = kk + 1;
                end
                tic;
                %[~, eomp_coeff, eomp_iopt, eomp_error,eomp_SSIM, eomp_hat]       = Enh_Convergence_OMP_pap7_(Dic, Vector,sparse,y_init,1, err_threshold(ii));
                [~, ~, ~,~, eomp_coeff,eomp_iopt,~, eomp_hat,eomp_error,eomp_SSIM] = Enh_OOMPex_PAP7(Vector, Dic, err_threshold(ii), sparse,y_init);
                [nocomb,E_SIZE_2,Approx_2,~] = runBRe(Vector,omp_error,eomp_coeff,eomp_iopt,Dic, Vector-eomp_hat(:,end), size(eomp_iopt,1), DicCoh,bits_per_index,bits_per_coeff,0,eomp_hat(:,end),y_init);
                %eomp_ci_matrix(eomp_iopt,kkk,ii,jj) = eomp_coeff;
                kkk = kkk + 1;
                EOMP_COEFF = [EOMP_COEFF;eomp_coeff];
                EOMP_IND = [EOMP_IND;eomp_iopt];

                %QC2 = imquantize(eomp_coeff,-5:0.0196:5,0:1:511);
                
                if mod(k,mod_threshold(jj)) == 0
                    y_init = eomp_hat(:,end).*0;
                else
                    y_init = eomp_hat(:,end);
                end
                T2 = toc;
                N2 = size(eomp_iopt,1);
                S2 = eomp_SSIM(end,1); zzzz = ssim_(Vector,Approx_2);zzzzz = ssim_(Vector,eomp_hat(:,end));
                E2 = eomp_error(end,1);
                BPP2 = (ceil(log2(size(omp_coeff,1))) + bits_per_coeff + (N2-1)*(bits_per_coeff + bits_per_index))/256;
                BPP2_BRe = E_SIZE_2/256;
                ALL_STHRE_MTHRE_NSEBT = [ALL_STHRE_MTHRE_NSEBT; err_threshold(ii), mod_threshold(jj), N1(i,j,ii),S1(i,j,ii),E1(i,j,ii),ori_bpp,BPP1(i,j,ii),BPP1_BRe(i,j,ii),T1(i,j,ii),N2,S2,E2,ori_bpp,BPP2,BPP2_BRe,T2];
                X_eOMP = [X_eOMP Vector];
                Xhat_eOMP = [Xhat_eOMP eomp_hat(:,end)];
                Xhat_eOMP_BRe = [Xhat_eOMP_BRe Approx_2];
                index = index + 1;
                k = k + 1;
            end
            if do_one_time == 0
                X_OMP_ = [X_OMP_; X_OMP];
                Xhat_OMP_ = [Xhat_OMP_; Xhat_OMP];
                Xhat_OMP_BRe_ = [Xhat_OMP_BRe_;Xhat_OMP_BRe];
            end
            X_eOMP_ = [X_eOMP_; X_eOMP];
            Xhat_eOMP_ = [Xhat_eOMP_; Xhat_eOMP];
            Xhat_eOMP_BRe_ = [Xhat_eOMP_BRe_; Xhat_eOMP_BRe];
        end
        
        [OMP_ALG1, OMP_ALG2] = AAA_getOptDenoised_OneStage_PAP7_(im2uint8(X_OMP_), im2uint8(Xhat_OMP_));
        [OMP_BRe_ALG1, OMP_BRe_ALG2] = AAA_getOptDenoised_OneStage_PAP7_(im2uint8(X_OMP_), im2uint8(Xhat_OMP_BRe_));
        [EOMP_ALG1, EOMP_ALG2] = AAA_getOptDenoised_OneStage_PAP7_(im2uint8(X_eOMP_), im2uint8(Xhat_eOMP_));
        [EOMP_BRe_ALG1, EOMP_BRe_ALG2] = AAA_getOptDenoised_OneStage_PAP7_(im2uint8(X_eOMP_), im2uint8(Xhat_eOMP_BRe_));
        
        SSIM_OMP = [ssim_(Xhat_OMP_,X_OMP_), ssim_(OMP_ALG1,X_OMP_),ssim_(OMP_ALG2,X_OMP_)];
        SSIM_OMP_BRe = [ssim_(Xhat_OMP_BRe_,X_OMP_), ssim_(OMP_BRe_ALG1,X_OMP_),ssim_(OMP_BRe_ALG2,X_OMP_)];
        SSIM_EOMP = [ssim_(Xhat_eOMP_,X_eOMP_), ssim_(EOMP_ALG1,X_eOMP_),ssim_(EOMP_ALG2,X_eOMP_)];
        SSIM_EOMP_BRe = [ssim_(Xhat_eOMP_BRe_,X_OMP_), ssim_(EOMP_BRe_ALG1,X_OMP_),ssim_(EOMP_BRe_ALG2,X_eOMP_)];
        PSNR_OMP = [psnr(im2uint8(Xhat_OMP_),im2uint8(X_OMP_)), psnr(im2uint8(OMP_ALG1),im2uint8(X_OMP_)),psnr(im2uint8(OMP_ALG2),im2uint8(X_OMP_))];
        PSNR_OMP_BRe = [psnr(im2uint8(Xhat_OMP_BRe_),im2uint8(X_OMP_)), psnr(im2uint8(OMP_BRe_ALG1),im2uint8(X_OMP_)),psnr(im2uint8(OMP_BRe_ALG2),im2uint8(X_OMP_))];
        PSNR_EOMP = [psnr(im2uint8(Xhat_eOMP_),im2uint8(X_eOMP_)), psnr(im2uint8(EOMP_ALG1),im2uint8(X_eOMP_)),psnr(im2uint8(EOMP_ALG2),im2uint8(X_eOMP_))];
        PSNR_EOMP_BRe = [psnr(im2uint8(Xhat_eOMP_BRe_),im2uint8(X_OMP_)), psnr(im2uint8(EOMP_BRe_ALG1),im2uint8(X_OMP_)),psnr(im2uint8(EOMP_BRe_ALG2),im2uint8(X_eOMP_))];
        NMSE_OMP = [immse(Xhat_OMP_,X_OMP_), immse(OMP_ALG1,X_OMP_),immse(OMP_ALG2,X_OMP_)];
        NMSE_OMP_BRe = [immse(Xhat_OMP_BRe_,X_OMP_), immse(OMP_BRe_ALG1,X_OMP_),immse(OMP_BRe_ALG2,X_OMP_)];
        NMSE_EOMP = [immse(Xhat_eOMP_,X_eOMP_), immse(EOMP_ALG1,X_eOMP_),immse(EOMP_ALG2,X_eOMP_)];
        NMSE_EOMP_BRe = [immse(Xhat_eOMP_BRe_,X_OMP_), immse(EOMP_BRe_ALG1,X_OMP_),immse(EOMP_BRe_ALG2,X_eOMP_)];
        


        [~,QOMP_COEFF_ind] = imquantize(OMP_COEFF,-5:0.0196:5,0:1:511);
        [~,QEOMP_COEFF_ind] = imquantize(EOMP_COEFF,-5:0.0196:5,0:1:511);
        [~, ~, QOMP_COEFF_bpp ] = get_prob_and_entropy(QOMP_COEFF_ind);
        [~, ~, QEOMP_COEFF_bpp ] = get_prob_and_entropy(QEOMP_COEFF_ind);
        [~, ~, OMP_IND_bpp ] = get_prob_and_entropy(OMP_IND);
        [~, ~, EOMP_IND_bpp ] = get_prob_and_entropy(EOMP_IND);
        
        BPP_ORI = ALL_STHRE_MTHRE_NSEBT(end,6);
        BPP_OMP = sum(ALL_STHRE_MTHRE_NSEBT(end-256+1:end,7))/256;
        BPP_OMP_BRe = sum(ALL_STHRE_MTHRE_NSEBT(end-256+1:end,8))/256- bits_per_index/256;
        BPP_OMP_ent = QOMP_COEFF_bpp + OMP_IND_bpp;
        BPP_OMP_BRe_ent = BPP_OMP_BRe * (BPP_OMP/BPP_OMP_ent)^(-1);

        %BPP_ORI = ALL_STHRE_MTHRE_NSEBT(end,6);
        BPP_EOMP = sum(ALL_STHRE_MTHRE_NSEBT(end-256+1:end,14))/256;
        BPP_EOMP_BRe = sum(ALL_STHRE_MTHRE_NSEBT(end-256+1:end,15))/256 - bits_per_index/256;% The index of the DC components is subtracted
        BPP_EOMP_ent = QEOMP_COEFF_bpp + EOMP_IND_bpp;
        BPP_EOMP_BRe_ent = BPP_EOMP_BRe * (BPP_EOMP/BPP_EOMP_ent)^(-1);

        OVERALL_BPP_RESULTS = [OVERALL_BPP_RESULTS; err_threshold(ii), mod_threshold(jj),...
            BPP_ORI,BPP_OMP,BPP_OMP_ent, BPP_OMP_BRe, BPP_OMP_BRe_ent, ...
            BPP_EOMP,BPP_EOMP_ent, BPP_EOMP_BRe, BPP_EOMP_BRe_ent];

        OVERALL_SSIM_RESULTS = [OVERALL_SSIM_RESULTS; err_threshold(ii), mod_threshold(jj),...
            SSIM_OMP,SSIM_OMP_BRe,SSIM_EOMP,SSIM_EOMP_BRe];
        
        OVERALL_PSNR_RESULTS = [OVERALL_PSNR_RESULTS; err_threshold(ii), mod_threshold(jj),...
            PSNR_OMP,PSNR_OMP_BRe,PSNR_EOMP,PSNR_EOMP_BRe];

        OVERALL_NMSE_RESULTS = [OVERALL_NMSE_RESULTS; err_threshold(ii), mod_threshold(jj),...
            NMSE_OMP,NMSE_OMP_BRe,NMSE_EOMP,NMSE_EOMP_BRe];


        do_one_time = 1;
        
        %OMP_C_P = get_prob(QOMP_COEFF_ind);
        %EOMP_C_P = get_prob(QEOMP_COEFF_ind);
        %OMP_I_P = get_prob(OMP_IND);
        %EOMP_I_P = get_prob(EOMP_IND);
        
        %OMP_COEFF_IND_LIST{ii,jj} = {QOMP_COEFF_ind OMP_IND};
        %EOMP_COEFF_IND_LIST{ii,jj} = {QEOMP_COEFF_ind EOMP_IND};
        
        check = 0;
    end
end
fprintf(strcat(img,' ... (Ours Alg. Reached %.2f%%)\n'),30*100/30);%2012
ALL_STHRE_MTHRE_NSEBT_mean = [];
for ii = 1:1:length(err_threshold)
    for jj = 1:1:length(mod_threshold)
        ind = find(ALL_STHRE_MTHRE_NSEBT(:,1) == err_threshold(ii) & ALL_STHRE_MTHRE_NSEBT(:,2) == mod_threshold(jj));
        ALL_STHRE_MTHRE_NSEBT_mean = [ALL_STHRE_MTHRE_NSEBT_mean; mean(ALL_STHRE_MTHRE_NSEBT(ind,:))];
    end
end
%%%%%
%%%%%OMP,OMP-BRe, -BM3D,-TVL1, EOMP,EOMP-BRe, -BM3D,-TVL1
%%%%%
%RESULTS.ALL_SSIM = 
%RESULTS.ALL_NMSE =
%RESULTS.ALL_SNR =
%RESULTS.ALL_bpp = 

RESULTS.ALL_STHRE_MTHRE_NSEBT = ALL_STHRE_MTHRE_NSEBT;
RESULTS.ALL_STHRE_MTHRE_NSEBT_mean = ALL_STHRE_MTHRE_NSEBT_mean;
RESULTS.OVERALL_SSIM_RESULTS = OVERALL_SSIM_RESULTS;
RESULTS.OVERALL_PSNR_RESULTS = OVERALL_PSNR_RESULTS;
RESULTS.OVERALL_NMSE_RESULTS = OVERALL_NMSE_RESULTS;
RESULTS.OVERALL_BPP_RESULTS = OVERALL_BPP_RESULTS;

RESULTS.JPEG_RES = JPEG_RES;
RESULTS.JPEG2000_RES = JPEG2000_RES;
RESULTS.ori_bpp = ori_bpp;
check = 0;
%RESULTS.eomp_ci_matrix = eomp_ci_matrix;
%RESULTS.omp_ci_matrix = omp_ci_matrix;
%RESULTS.ORI_IMG_double = img_double;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SS_den, NM_den, imgs] = buildImages(yi_omp,yi_iomp,yi_nomp, r, c,img_pixels,img_double)
img_omp = zeros(r*8, c*8,12);
img_iomp = zeros(r*8, c*8,12);
img_nomp = zeros(r*8, c*8,12);


for sp = 1:1:12
    fprintf('Processing sp %.2f >>> %.2f %%',sp);
    x_c = [];
    y_c = [];
    z_c = [];
    for i = 1:1:r
        x_r = [];
        y_r = [];
        z_r = [];
        for j=1:1:c
            
            
            v = yi_omp(1:end,sp,i,j);
            omp_block = vectorToBatch(v,8);
            v = yi_iomp(1:end,sp,i,j);
            iomp_block = vectorToBatch(v,8);
            v = yi_nomp(1:end,sp,i,j);
            nomp_block = vectorToBatch(v,8);
            %end
            x_r = [x_r omp_block];
            y_r = [y_r iomp_block];
            z_r = [z_r nomp_block];
            stop_here = 0;
            
            
            
        end
        x_c = [x_c; x_r];
        y_c = [y_c; y_r];
        z_c = [z_c; z_r];
        stop_here = 0;
    end
    img_omp(:,:,sp) = x_c;
    img_iomp(:,:,sp) = y_c;
    img_nomp(:,:,sp) = z_c;
    
    img_omp_im(:,:,sp) = im2uint8(x_c);
    img_iomp_im(:,:,sp) = im2uint8(y_c);
    img_nomp_im(:,:,sp) = im2uint8(z_c);
    
    s1 = ssim(im2uint8(img_double(:,:,1)),im2uint8(x_c));
    s2 = ssim(im2uint8(img_double(:,:,1)),im2uint8(y_c));
    s3 = ssim(im2uint8(img_double(:,:,1)),im2uint8(z_c));
    
    %figure
    subplot(1,3,1)
    imshow(im2uint8(x_c));
    subplot(1,3,2)
    imshow(im2uint8(y_c));
    subplot(1,3,3)
    imshow(im2uint8(z_c));
    cc = 0;
    clc;
    
end
SS_den.omp = SS_omp_den;
SS_den.iomp = SS_iomp_den;
SS_den.nomp = SS_nomp_den;
NM_den.omp = NM_omp_den;
NM_den.iomp = NM_iomp_den;
NM_den.nomp = NM_nomp_den;

imgs.omp = img_omp_im;
imgs.iomp = img_iomp_im;
imgs.nomp = img_nomp_im;
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
function patch = vec2Patch(vec)
if vec(1)>=1
    patch(1:8,1) = vec(1:8)./255;
    patch(1:8,2) = vec(9:16)./255;
    patch(1:8,3) = vec(17:24)./255;
    patch(1:8,4) = vec(25:32)./255;
    patch(1:8,5) = vec(33:40)./255;
    patch(1:8,6) = vec(41:48)./255;
    patch(1:8,7) = vec(49:56)./255;
    patch(1:8,8) = vec(57:64)./255;
else
    patch(1:8,1) = vec(1:8);
    patch(1:8,2) = vec(9:16);
    patch(1:8,3) = vec(17:24);
    patch(1:8,4) = vec(25:32);
    patch(1:8,5) = vec(33:40);
    patch(1:8,6) = vec(41:48);
    patch(1:8,7) = vec(49:56);
    patch(1:8,8) = vec(57:64);
end
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

% -------------------------------------------------------------------------
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

function       [NoCombs_2,SIZE_2,Approx_,TIME_2] = runBRe(X,EOMP,COEFF,IOPT,Dic, R, Sparsity_, DicCoh,bitsPerIndex,bitsPerCoeff,TIME_2,YYY_2,y_init)

[TIME2, NoCombs_2,Approx_]=ERR_based_BRe_modified(X,EOMP,COEFF,IOPT,Dic, R, Sparsity_, DicCoh, y_init);
%YYY_2 = [YYY_2;X-ErrMat2(:,end)];
SIZE_2 =  Sparsity_ * bitsPerIndex + (Sparsity_ - NoCombs_2)*bitsPerCoeff;
%CR_BRe = (SIZE_1-SIZE_2)*100/SIZE_1;
TIME_2 = TIME_2 + TIME2;
end


function [time1, NoComb,Approx_] = ERR_based_BRe_modified(Y, errFS, coeffFS, indicesFS, Dic, rFS, Spa, DicCoh, y_init )
%build k*k matrix for 0.5*(Ci-Cj)^2
%for +ve and -ve coeffs 
tic;

if isempty(y_init)
    y_init = zeros(size(Y,1),1);
end
CoFs = coeffFS;
InFs = indicesFS;

%fprintf('Alg. %s\n','BaSp');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TWOS = [];
for iii = 1:1:Spa-1
    for jjj = iii+1:1:Spa
        %for kkk = jjj+1:1:Spa
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Compute Coefficient
            w_ = (CoFs(iii)+CoFs(jjj))/2;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Compute Error
            
            error_ = 0.5*(CoFs(iii)-CoFs(jjj))^2*(1-DicCoh(InFs(iii),InFs(jjj))) ;



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            TWOS = [TWOS; InFs(iii) InFs(jjj) w_ error_];
        %end
    end
end
SELECTED_COMBS = [];
for sel = 1:1:floor(Spa/2)
    [C2,I2] = min(TWOS(:,4));
    selected = TWOS(I2,:);
    SELECTED_COMBS = [SELECTED_COMBS; TWOS(I2,:)];
    toRemove = [];
    for x = 1:1:size(TWOS,1)
        if(TWOS(x,1) == selected(1) || TWOS(x,1) == selected(2))
            toRemove = [toRemove, x];
        end

    end
    TWOS(toRemove,:) = [];

    toRemove = [];
    for x = 1:1:size(TWOS,1)
        if(TWOS(x,2) == selected(1) || TWOS(x,2) == selected(2))
            toRemove = [toRemove, x];
        end

    end
    TWOS(toRemove,:) = [];

end



%Y, yFS, errFS, coeffFS, indicesFS, Dic, rFS, Spa, DicCoh, Code
COEFFS = zeros(size(Dic,2),1);
COEFFS(indicesFS,1)=coeffFS;
NoComb = 0;
ErrMat1 = Y - Dic *COEFFS- y_init;
dd = norm(ErrMat1);
COEFFS_FIN = COEFFS;
try
for i=1:1:size(SELECTED_COMBS,1)
    COEFFS([SELECTED_COMBS(i,1),SELECTED_COMBS(i,2)],1) = SELECTED_COMBS(i,3);
     %ee = Y - Dic *COEFFS;
     %dddd = norm(ee);
     %check = 0;
    if immse(Y , Dic *COEFFS + y_init) <= errFS(end-1)
        NoComb = NoComb + 1;
        ErrMat1(:,i) = Y - Dic *COEFFS - y_init;
        COEFFS_FIN = COEFFS;
    else

        break;
    end
    

end
catch
end
        Approx_ = Dic *COEFFS_FIN + y_init;
        %xxcx = ssim_(Y, Dic *COEFFS_FIN + y_init);
        %xxcxx = ssim_(Y,y_init);

cc = 0;
time1 = toc;

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


function [ D, Di, beta,betaa, CC, Di_,ErrMat,APP,Err_,oomp_ssim] = OOMPex_PAP7( f, D, tol, No, ind )
% OOMP  Optimized Orthogonal Matching Pursuit
%
% It creates an atomic decomposition of a signal using OOMP method [1]. You can choose a
% tolerance, the number of atoms to take in or an initial subspace to influence the OOMP
% algorithm. Non-selected atoms subtracted by their component in the chosen space are also
% available.
%
% Usage:    [ Dnew, beta, Di ] = OOMP( f, D, tol );
%           [ Dnew, beta, Di, c, Q ] = OOMP( f, D, tol, No, ind );
%
% Inputs:
%   f       signal to be represented
%   D       dictionary of atoms
%   tol     desired distance between f and its approximation the routine will stop if
%           norm(f'-Dsub*(f*beta)')*sqrt(delta)<tol where delta=1/L, L is number of points
%           in a sample
%   No      (optional) maximal number of atoms to choose, if the number of chosen atoms
%           equals to No, routine will stop (default No=size(D,2))
%   ind     (optional) indices determining  the initial subspace,
%
% Outputs:
%   D       the dictionary D rearranged according to the selection process D(:,1:k)
%           contains the atoms chosen into the atomic decomposition
%   beta    'k' biorthogonal functions corresponding to new D(:,1:k)
%   Di      indices of atoms in new D written w.r.t the original D
%   c       'k' coefficients of the atomic decomposition
%   Q       Q(:,1:k) contains orthonormal functions spanning new D(:,1:k), Q(:,k+1:N)
%           contains new D(:,k+1:N) subtracted by the projection onto the space generated
%           by  Q(:,1:k)  (resp. D(:,1:k))
%
% References:
%   [1] L. Rebollo-Neira and D. Lowe, "Optimized Orthogonal Matching Pursuit Approach",
%   IEEE Signal Processing Letters, Vol(9,4), 137-140, (2002).
%   For the current implementation:
%   [2] M. Andrle and L. Rebollo-Neira, "A swapping-based refinement of orthogonal
%   matching pursuit strategies", Signal Processing, Vol 86, No 3, pp. 480-495, (2006).
%
% See also OMP Swapping OOMPKSwaps OMPKSwaps OOMPFinalRefi BOOMP OBOMP

% See   http://www.ncrg.aston.ac.uk/Projects/HNLApprox/
%       http://www.nonlinear-approx.info/
name='OOMP';  %name of routine
% get inputs and setup parameters
Err=[zeros(No,1)];
indices = [zeros(No,1)];
[L,N]=size(D);
beta=[];
Di=1:N;     % index vector of full-dictionary
Q=D;
%delta=1/L;  %uncomment in the case of analytical norm
delta=1;   %uncomment in the case of discrete norm
if nargin<5 ind=[];end
if nargin<4 No=N;end
if nargin<3 tol=0.01;end;

numind=length(ind);

%atoms having smaller norm than tol1 are supposed be zero ones
tol1=1e-7; %1e-5
%threshold for coefficients
tol2=1e-20;%tol2=1e-10;   %0.0001  %1e-5

%tic;
%fprintf('\n%s is running for tol=%g, tol1=%g and tol2=%g.\n',name,tol,tol1,tol2);
%fprintf('tol  -> required precision, distance ||f-f_approx||\n')
%fprintf('tol1 -> atoms having smaller norm are supposed to be zero ones.\n');
%fprintf('tol2 -> algorithm stops if max(|<f,q>|/||q||)<= tol2.\n');
%===============================
% Main algorithm: at kth iteration
%===============================
H=min(No,N); %maximal number of function in sub-dictionary
for k=1:H
    %finding of maximal coefficient
    c=zeros(1,N);cc=zeros(1,N);
    if k<=numind
        [test,q]=ismember(ind(k),Di);

        if test~=1; error('Demanded index %d is out of dictionary',ind(k));end
        c(k)=norm(Q(:,q));
        if c(k)<tol1
            error('Demanded atom with index %d is dependent on the others.',ind(k));
        end;
    else
        for p=k:N, c(p)=norm(Q(:,p));
            if c(p)>tol1; cc(p)=abs(f'*Q(:,p))/c(p);end;
        end
        [max_c,q]=max(cc);
        indices(k,1) = q;
        %stopping criterion (coefficient)
        if max_c<tol2 ;
            k=k-1;
            %fprintf('%s stopped, max(|<f,q>|/||q||)<= tol2=%g.\n',name,tol2);
            break;
        end
    end
    if q~=k;
        Q(:,[k q])=Q(:,[q k]); % swap k-th and q-th columns
        D(:,[k q])=D(:,[q k]);
        Di([k q])=Di([q k]);
    end
    %re-orthogonalization of Q(:,k)  w.r.t Q(:,1),..., Q(:,k-1)
    if k>1;
        for p=1:k-1;
            Q(:,k)=Q(:,k)-(Q(:,p)'*Q(:,k))*Q(:,p);
        end
    end
    nork=norm(Q(:,k));
    Q(:,k)=Q(:,k)/nork; %normalization
    % compute biorthogonal functions beta from 1 to k-1
    if k>1;
        beta=beta-Q(:,k)*(D(:,k)'*beta)/nork;
    end
    beta(:,k)=Q(:,k)/nork;          % kth biorthogonal function
    betaa{k} = beta;
    %orthogonalization of Q(:,k+1:n) wrt Q(:,k)
    if k==N, break; end
    for p=k+1:N;
        % orthogonalization of Q(:,p) wrt Q(:,k)
        Q(:,p)=Q(:,p)-(Q(:,k)'*Q(:,p))*Q(:,k);
    end
    %stopping criterion (distance)
    app = D(:,1:k)*(f'*beta)';
    APP(:,k)=app;
    xxxxx = f-app;
    %Err_(k,1) = norm(xxxxx);
    Err_(k,1) = immse(f,app);
    ErrMat(:,k)=xxxxx;
    if (tol~= 0) && (norm(xxxxx)*sqrt(delta)*100/norm(f) < tol); break;end;
    %SNRx(k,1)=snr(f,xxxxx);
    %Err(k,1)=norm(xxxxx)*sqrt(delta);
    vvvv=6;
    %fprintf('Alg. %s, ite. %i\n','OOMP',k);


    %try
    %[PesqMos] = comp_pesq_GLOBE(f,16000,app,16000);
    %catch
        %PesqMos = [1 1];
    %end
 

end

c=f'*beta;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normf=norm(f-D(:,1:k)*c')*sqrt(delta);
%ROOMP=f-D(:,1:k)*c';
oomp_ssim = ssim_(D(:,1:k)*c',f);
%fprintf('From %d atoms in this dictionary %s has chosen %d atoms.\n',N,name,k);
%fprintf('\n%s lasted %g seconds\n',name,toc);
%fprintf('\nNorm ||f-f_approx|| is %g. \n',normf);
% error testing
% orthogonality
ErrorTest(Q(:,1:k));
% biorthogonality
ErrorTest(D(:,1:k),beta);
%time = toc;
CC = c';
Di_ = Di(1:k')';
%Last modification Laura REBOLLO-NEIRA (2009) (former OOMPF)
%
%Copyright (C) 2006 Miroslav ANDRLE and Laura REBOLLO-NEIRA
%
%This program is free software; you can redistribute it and/or modify it under the terms
%of the GNU General Public License as published by the Free Software Foundation; either
%version 2 of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
%without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%See the GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License along with this program;
%if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
%Boston, MA  02110-1301, USA.
cc = 0;
end

function [ D, Di, beta,betaa, CC, Di_,ErrMat,APP,Err_,oomp_ssim] = Enh_OOMPex_PAP7( f, D, tol, No, y_init )
% OOMP  Optimized Orthogonal Matching Pursuit
%
% It creates an atomic decomposition of a signal using OOMP method [1]. You can choose a
% tolerance, the number of atoms to take in or an initial subspace to influence the OOMP
% algorithm. Non-selected atoms subtracted by their component in the chosen space are also
% available.
%
% Usage:    [ Dnew, beta, Di ] = OOMP( f, D, tol );
%           [ Dnew, beta, Di, c, Q ] = OOMP( f, D, tol, No, ind );
%
% Inputs:
%   f       signal to be represented
%   D       dictionary of atoms
%   tol     desired distance between f and its approximation the routine will stop if
%           norm(f'-Dsub*(f*beta)')*sqrt(delta)<tol where delta=1/L, L is number of points
%           in a sample
%   No      (optional) maximal number of atoms to choose, if the number of chosen atoms
%           equals to No, routine will stop (default No=size(D,2))
%   ind     (optional) indices determining  the initial subspace,
%
% Outputs:
%   D       the dictionary D rearranged according to the selection process D(:,1:k)
%           contains the atoms chosen into the atomic decomposition
%   beta    'k' biorthogonal functions corresponding to new D(:,1:k)
%   Di      indices of atoms in new D written w.r.t the original D
%   c       'k' coefficients of the atomic decomposition
%   Q       Q(:,1:k) contains orthonormal functions spanning new D(:,1:k), Q(:,k+1:N)
%           contains new D(:,k+1:N) subtracted by the projection onto the space generated
%           by  Q(:,1:k)  (resp. D(:,1:k))
%
% References:
%   [1] L. Rebollo-Neira and D. Lowe, "Optimized Orthogonal Matching Pursuit Approach",
%   IEEE Signal Processing Letters, Vol(9,4), 137-140, (2002).
%   For the current implementation:
%   [2] M. Andrle and L. Rebollo-Neira, "A swapping-based refinement of orthogonal
%   matching pursuit strategies", Signal Processing, Vol 86, No 3, pp. 480-495, (2006).
%
% See also OMP Swapping OOMPKSwaps OMPKSwaps OOMPFinalRefi BOOMP OBOMP

% See   http://www.ncrg.aston.ac.uk/Projects/HNLApprox/
%       http://www.nonlinear-approx.info/
name='OOMP';  %name of routine
% get inputs and setup parameters
Err=[zeros(No,1)];
indices = [zeros(No,1)];
[L,N]=size(D);
beta=[];
Di=1:N;     % index vector of full-dictionary
Q=D;
%delta=1/L;  %uncomment in the case of analytical norm
delta=1;   %uncomment in the case of discrete norm
%if nargin<5 ind=[];end
%if nargin<4 No=N;end
%if nargin<3 tol=0.01;end;
ind = [];

numind=length(ind);
%No=N;
%No = 128;
%atoms having smaller norm than tol1 are supposed be zero ones
tol1=1e-7; %1e-5
%threshold for coefficients
tol2=1e-20;%tol2=1e-10;   %0.0001  %1e-5
ff = f - y_init;
%f = ff;
%tic;
%fprintf('\n%s is running for tol=%g, tol1=%g and tol2=%g.\n',name,tol,tol1,tol2);
%fprintf('tol  -> required precision, distance ||f-f_approx||\n')
%fprintf('tol1 -> atoms having smaller norm are supposed to be zero ones.\n');
%fprintf('tol2 -> algorithm stops if max(|<f,q>|/||q||)<= tol2.\n');
%===============================
% Main algorithm: at kth iteration
%===============================
H=min(No,N); %maximal number of function in sub-dictionary
for k=1:H
    %finding of maximal coefficient
    c=zeros(1,N);cc=zeros(1,N);
    if k<=numind
        [test,q]=ismember(ind(k),Di);

        if test~=1; error('Demanded index %d is out of dictionary',ind(k));end
        c(k)=norm(Q(:,q));
        if c(k)<tol1
            error('Demanded atom with index %d is dependent on the others.',ind(k));
        end
    else
        for p=k:N, c(p)=norm(Q(:,p));
            if c(p)>tol1; cc(p)=abs(ff'*Q(:,p))/c(p);end
        end
        [max_c,q]=max(cc);
        indices(k,1) = q;
        %stopping criterion (coefficient)
        if max_c<tol2 
            k=k-1;
            %fprintf('%s stopped, max(|<f,q>|/||q||)<= tol2=%g.\n',name,tol2);
            break;
        end
    end
    if q~=k
        Q(:,[k q])=Q(:,[q k]); % swap k-th and q-th columns
        D(:,[k q])=D(:,[q k]);
        Di([k q])=Di([q k]);
    end
    %re-orthogonalization of Q(:,k)  w.r.t Q(:,1),..., Q(:,k-1)
    if k>1
        for p=1:k-1
            Q(:,k)=Q(:,k)-(Q(:,p)'*Q(:,k))*Q(:,p);
        end
    end
    nork=norm(Q(:,k));
    Q(:,k)=Q(:,k)/nork; %normalization
    % compute biorthogonal functions beta from 1 to k-1
    if k>1
        beta=beta-Q(:,k)*(D(:,k)'*beta)/nork;
    end
    beta(:,k)=Q(:,k)/nork;          % kth biorthogonal function
    betaa{k} = beta;
    %orthogonalization of Q(:,k+1:n) wrt Q(:,k)
    if k==N, break; end
    for p=k+1:N
        % orthogonalization of Q(:,p) wrt Q(:,k)
        Q(:,p)=Q(:,p)-(Q(:,k)'*Q(:,p))*Q(:,k);
    end
    %stopping criterion (distance)
    app = D(:,1:k)*(ff'*beta)';
    APP(:,k)=app+y_init;
    xxxxx = ff-app;
    %Err_(k,1) = norm(xxxxx);
    Err_(k,1) = immse(f,app+y_init);
    ErrMat(:,k)=xxxxx;
    %temp = y - (e_hat +y_init) ;
    %curERR(count,1) = norm(temp);
    %err_ = curERR(count,1)*100/normy;
    zzz  = (norm(f - (app + y_init))*100/norm(f));
    if (tol~= 0) && (zzz < tol) 
        break;
    end
    %SNRx(k,1)=snr(f,xxxxx);
    %Err(k,1)=norm(xxxxx)*sqrt(delta);
    vvvv=6;
    %fprintf('Alg. %s, ite. %i\n','OOMP',k);


    %try
    %[PesqMos] = comp_pesq_GLOBE(f,16000,app,16000);
    %catch
        %PesqMos = [1 1];
    %end
 

end

c=ff'*beta;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normf=norm(f-D(:,1:k)*c')*sqrt(delta);
%ROOMP=f-D(:,1:k)*c';
oomp_ssim = ssim_(f,APP(:,end));
%fprintf('From %d atoms in this dictionary %s has chosen %d atoms.\n',N,name,k);
%fprintf('\n%s lasted %g seconds\n',name,toc);
%fprintf('\nNorm ||f-f_approx|| is %g. \n',normf);
% error testing
% orthogonality
ErrorTest(Q(:,1:k));
% biorthogonality
ErrorTest(D(:,1:k),beta);
%time = toc;
CC = c';
Di_ = Di(1:k')';
%Last modification Laura REBOLLO-NEIRA (2009) (former OOMPF)
%
%Copyright (C) 2006 Miroslav ANDRLE and Laura REBOLLO-NEIRA
%
%This program is free software; you can redistribute it and/or modify it under the terms
%of the GNU General Public License as published by the Free Software Foundation; either
%version 2 of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
%without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%See the GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License along with this program;
%if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
%Boston, MA  02110-1301, USA.
cc = 0;
end


function f = ErrorTest( D, beta )
% ErrorTest tests orthogonality of a sequence or biorthogonality of two sequences.
% 
% Usage:    errortest( D, beta );
%           errortest( Q );
%
% Inputs:
%   D       sequence of vectors 
%   beta    (optional) biorthogonal sequence  
%
% Output:
%   f       orthogonality of D or biorthogonality D w.r.t beta

% See   http://www.ncrg.aston.ac.uk/Projects/HNLApprox/
%       http://www.nonlinear-approx.info/

name='ErrorTest';  
  
if nargin==1   
  % orthogonality
  PP=D'*D;f=norm(eye(size(PP))-PP);
  %fprintf('Orthogonality:   %g \n',f);
elseif nargin==2
  % biorthogonality
  PP=beta'*D;f=norm(eye(size(PP))-PP);
  %fprintf('Biorthogonality: %g \n',f);
else fprintf('%s: Wrong number of arguments',name);  
end  


%Copyright (C) 2006 Miroslav ANDRLE and Laura REBOLLO-NEIRA
%
%This program is free software; you can redistribute it and/or modify it under the terms 
%of the GNU General Public License as published by the Free Software Foundation; either 
%version 2 of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
%without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%See the GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License along with this program;
%if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
%Boston, MA  02110-1301, USA.
end