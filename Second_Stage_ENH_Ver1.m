
function [RESULTS] = Second_Stage_ENH_Ver1(img,batchSize,sparse_,dicType) %#ok<*INUSD>
%BRE_IMG Summary of this function goes here
%   Detailed explanation goes here
global sparse;
sparse = sparse_;
img_pixels = imread(img);
temp_comp_ssim = zeros(12,3);
iter = 0;
vecLength = batchSize*batchSize;
img_double = im2double(img_pixels);
Dic = getDic(dicType,vecLength);
[Batches, batchPerRow,batchPerCol] = getBatches(img_double,batchSize);
global noisyIMAGES;
global sOMP_noisyIMAGES;
noisyIMAGES = zeros(batchSize*batchPerCol,batchSize*batchPerRow,sparse+1);%to build the overall OMP-based images, first image is the origional one
sOMP_noisyIMAGES = zeros(batchSize*batchPerCol,batchSize*batchPerRow,sparse+1);%to build the overall sOMP-based images, first image is the origional one

index = 1;
COLLECTED_SSIM = zeros(sparse,6);
COLLECTED_NMSE = zeros(sparse,6);
k=1;
algs = 6;
ROS_m = zeros(sparse,algs); 
RON_m = zeros(sparse,algs); 


xxx =[];
yyy = [];
zzz = [];
sig = [];
Delay_per_ite = [];
Delay_per_ite_locbe = [];
for i=1:1:floor(batchPerRow)
    for j=1:1:floor(batchPerCol)
        Vector = batchToVector(Batches(i,j,:,:),batchSize);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%(OMP, sOMP)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Ignored Output
        %RESULTS_nOMP(k) RESULTS_iOMP(k) RESULTS_sOMP(k) RESULTS_LocBE(k) RESULTS_FBP(k)
      
        [~, coeff, iopt, OMP_curERR,omp_ssim, y_i_omp] = SSIM_OMP_new(Dic, Vector,sparse);
        [~, iomp_curERR, iomp_ssim, y_iOMP]                      = SSIM_iOMP_new(Dic, Vector,sparse);   
        iOMP_ssim = iomp_ssim';
        [miomp_ssim, miomp_curERR, y_i_miomp]                      = getMIOMP(omp_ssim, iomp_ssim,OMP_curERR, iomp_curERR,y_i_omp,y_iOMP, sparse); 
        MiOMP_ssim = miomp_ssim';
        [~, somp_curERR, somp_ssim, y_sOMP]                = SSIM_selection_OMP_new(Dic, Vector,sparse);             
        [~, coeff, iopt, OMP_curERR,...
        nOMP_curERR,omp_ssim,nomp_ssim, y_i_omp, y_i_nomp,ssimoOPT_time]       = SSIM_OMP_nOMP_new(Dic, Vector,sparse);
        nOMP_ssim = nomp_ssim';
        %Update Noisy Images
        %noisyIMAGES = updateNoisyIMG(Vector,y_i_omp,batchSize,batchPerCol,batchPerRow,i,j,sparse,noisyIMAGES);

        %[~, somp_curERR, somp_ssim, y_sOMP]                = SSIM_selection_OMP_new(Dic, Vector,sparse);             
        %sOMP_noisyIMAGES = updateNoisyIMG(Vector,y_sOMP,batchSize,batchPerCol,batchPerRow,i,j,sparse,sOMP_noisyIMAGES);
        stop = 0;
        %[~, LocBE_curERR,  LocBE_ssim, y_i_LocBE]                = LocBE_new(Dic, Vector,coeff,iopt,omp_ssim,y_i_omp);
        %sOMP_noisyIMAGES = updateNoisyIMG(Vector,y_i_LocBE,batchSize,batchPerCol,batchPerRow,i,j,sparse,sOMP_noisyIMAGES);

            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
         
        %for kk =1:1:sparse
            %[SS_omp_den(kk,:) NM_omp_den(kk,:) bm3d_runtime1] = getOptDenoised_OneStage(Vector, y_i_omp(:,kk),omp_ssim(kk),OMP_curERR(kk));
        %end

         %ROS(:,:,k) =  SS_omp_den;
         %RON(:,:,k) =  NM_omp_den;

         
         
         %ROS_m  =  ROS_m + ROS(:,:,k);
         %RON_m  =  RON_m + RON(:,:,k);

         
         %COLLECTED_SSIM = COLLECTED_SSIM + [omp_ssim'];
         %COLLECTED_NMSE = COLLECTED_NMSE + [OMP_curERR];
         %COLLECTED_SSIM_avg = COLLECTED_SSIM./k;
         %COLLECTED_NMSE_avg = COLLECTED_NMSE./k;
         %Check = 0;

          clc;
          fprintf('Processing %s >>> %.2f %%',img,index*100/(batchPerRow*batchPerCol));

          index = index + 1;  
          
          
          k = k + 1;
     end
end
        %for kk =1:1:sparse
            %[SS_omp_den(kk,:) NM_omp_den(kk,:) bm3d_runtime1] = getOptDenoised_OneStage(Vector, y_i_omp(:,kk),omp_ssim(kk),OMP_curERR(kk));
[SS_omp_den(kk,:) NM_omp_den(kk,:) bm3d_runtime1] = getDenoisedIMAGE(noisyIMAGES,sOMP_noisyIMAGES,sparse);
        %end
ROS_m  =  ROS_m./k;
RON_m  =  RON_m./k;
check  = 0;
RESULTS.ROS  =  ROS;
RESULTS.RON  =  RON;

RESULTS.ROS_m  =  ROS_m;
RESULTS.RON_m  =  RON_m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% we do the for_back processing for the main algs%%%%%%%%
%%%%%%%%% OMP, iOMP, FBP, MiOMP, nOMP, LocBE, sOMP%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[SS_den, NM_den, imgs] = buildImages(y_i_omp_,y_i_iomp_,y_i_nomp_,floor(batchPerRow),floor(batchPerCol),img_pixels,img_double);

RESULTS.FOR_RESULTS = forProcessing(RESULTS, sparse,img);
RESULTS.BAC_RESULTS = bacProcessing(RESULTS, sparse,img);
         
check = 0;

end

function   noisyIMAGES = updateNoisyIMG(Vector,y_i_omp,batchSize,batchPerCol,batchPerRow,i,j,sparse,noisyIMAGES)
batchRowIndex = i;
batchColIndex = j;
Upper_Left_x = batchSize*(batchColIndex-1)+1;
Upper_Left_y = batchSize*(batchRowIndex-1)+1;
for k = 1:1:sparse+1
    element = 1;
    for row = Upper_Left_y:1:Upper_Left_y + batchSize - 1
        for col = Upper_Left_x:1:Upper_Left_x + batchSize - 1
            if k == 1
                noisyIMAGES(row,col,k) = Vector(element,k);
            else
                noisyIMAGES(row,col,k) = y_i_omp(element,k-1);
            end
            
            element = element + 1;
        end
    end
end
xx = noisyIMAGES(:,:,1);
stop = 0;
end

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
 
            %omp_block = zeros(8,8);
            %iomp_block = zeros(8,8);
            %nomp_block = zeros(8,8);    
            %for ii= 1:1:8
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
   %[SS_omp_den(sp,:) NM_omp_den(sp,:)] = getOptDenoised_OneStage_nOMP(img_pixels(:,:,1),img_double(:,:,1), x_c, im2uint8(x_c));
   %[SS_iomp_den(sp,:) NM_iomp_den(sp,:)] = getOptDenoised_OneStage_nOMP(img_pixels(:,:,1),img_double(:,:,1), y_c, im2uint8(y_c));
   %[SS_nomp_den(sp,:) NM_nomp_den(sp,:)] = getOptDenoised_OneStage_nOMP(img_pixels(:,:,1),img_double(:,:,1), z_c, im2uint8(z_c));
   
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










