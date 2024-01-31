function [RESULTS] = Second_BE_Final_SparseALGs(img,batchSize,sparse_,dicType) %#ok<*INUSD>
%BRE_IMG Summary of this function goes here
%   Detailed explanation goes here
global sparse;
sparse = sparse_;
img_pixels = imread(img);
iter = 0;
vecLength = batchSize*batchSize;
img_double = im2double(img_pixels);
Dic = getDic(dicType,vecLength);

[Batches, batchPerRow,batchPerCol] = getBatches(img_double,batchSize);
index = 1;
COLLECTED_SSIM = zeros(sparse,4);
COLLECTED_NMSE = zeros(sparse,4);
k=1;
algs = 6;
ROS_m = zeros(sparse,algs); 
RON_m = zeros(sparse,algs); 
RiS_m = zeros(sparse,1); 
RiN_m = zeros(sparse,1); 
%RFS_m = zeros(sparse,1); 
%RFN_m = zeros(sparse,1); 
RmS_m = zeros(sparse,algs);  
RmN_m = zeros(sparse,algs);  
RnS_m = zeros(sparse,algs);  
RnN_m = zeros(sparse,algs);  
RsS_m = zeros(sparse,1);  
RsN_m = zeros(sparse,1);  
%RLS_m = zeros(sparse,1);  
%RLN_m = zeros(sparse,1); 
xxx =[];
yyy = [];
zzz = [];
sig = [];
Delay_per_ite = [];
Delay_per_ite_locbe = [];
SSIM_ALL = zeros(sparse,12);
ERROR_ALL = zeros(sparse,12);
warning('off','all');
y_init = zeros(batchSize^2,1);
supp_init =[];
for i=1:1:floor(batchPerRow)
    for j=1:1:floor(batchPerCol)
        Vector = batchToVector(Batches(i,j,:,:),batchSize);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(OMP(1990) - StOMP(2006)? - SWOMP(2009)  LASSO(2009) - COSAMP(2008) - GSM(2020)-
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OGA() - OMPvar(2015) - SAEN(2018) - GOMP(2010) -  - cpwwen(2018) - clarswlasso(2018) )%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(strcat('Processing ' , img , ' by >>>' , ' (OMP)\n'));
        [~, coeff, iopt, OMP_ERROR,OMP_SSIM, y_i_omp]       = SSIM_OMP_new(Dic, Vector,sparse);
        
        [~, coeff, iopt_, OMP_ERROR_2,OMP_SSIM_2, y_i_omp_n]       = Enh_Convergence_OMP(Dic, Vector,sparse,y_init,supp_init);
        y_init = y_i_omp_n(:,end);
        supp_init = iopt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        fprintf(strcat('Processing ' , img , ' by >>>' , ' (iOMP)'));%2012
        [~, iOMP_ERROR, iOMP_SSIM, y_iOMP]                      = SSIM_iOMP_new(Dic, Vector,sparse);  
        RES = find(iOMP_SSIM(:)>OMP_SSIM(:));
        RESU = find(iOMP_ERROR(:)<OMP_ERROR(:));        
        fprintf(strcat(', (%d%%, %d%%)\n'),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
        fprintf(strcat('Processing ' , img , ' by >>>' , ' (MiOMP)'));%2021
        [MiOMP_SSIM, MiOMP_ERROR, y_i_miomp]                      = getMIOMP(OMP_SSIM, iOMP_SSIM,OMP_ERROR, iOMP_ERROR,y_i_omp,y_iOMP, sparse); 
        RES = find(MiOMP_SSIM(:)>OMP_SSIM(:));
        RESU = find(MiOMP_ERROR(:)<OMP_ERROR(:));                
        fprintf(strcat(', (%d%%, %d%%)\n'),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
        fprintf(strcat('Processing ' , img , ' by >>>' , ' (nOMP)'));%2022
        [~, coeff, iopt, OMP_curERR, nOMP_ERROR,omp_ssim,nOMP_SSIM, y_i_omp, y_i_nomp,ssimoOPT_time]       = SSIM_OMP_nOMP_new(Dic, Vector,sparse);
        RES = find(nOMP_SSIM(:)>OMP_SSIM(:));
        RESU = find(nOMP_ERROR(:)<OMP_ERROR(:));                        
        fprintf(strcat(', (%d%%, %d%%)\n'),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        fprintf(strcat('Processing ' , img , ' by >>>' , ' (LocBE)'));
        [~, LocBE_ERROR,  LocBE_SSIM, y_i_LocBE]                = LocBE_new(Dic, Vector,coeff,iopt,omp_ssim,y_i_omp);
        RES = find(LocBE_SSIM(:)>OMP_SSIM(:));
        RESU = find(LocBE_ERROR(:)<OMP_ERROR(:));                        
        fprintf(strcat(', (%d%%, %d%%)\n'),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        fprintf(strcat('Processing ' , img , ' by >>>' , ' (SAMP)'));       
        [xr, iter_num,SAMP_SSIM,SAMP_ERROR] = SAMP_ano(Vector, Dic, 1, sparse,0);%2008
        RES = find(SAMP_SSIM(:)>OMP_SSIM(:));
        RESU = find(SAMP_ERROR(:)<OMP_ERROR(:));                                
        fprintf(strcat(', (%d%%, %d%%)\n'),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        fprintf(strcat('Processing ' , img , ' by >>>' , ' (COSAMP)'));
        [Sest,COSAMP_ERROR,COSAMP_SSIM] = cosamp_ano(Dic,Vector,sparse,0,20);%2009
        RES = find(COSAMP_SSIM(:)>OMP_SSIM(:));
        RESU = find(COSAMP_ERROR(:)<OMP_ERROR(:));                                        
        fprintf(strcat(', (%d%%, %d%%)\n'),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        
       fprintf(strcat('Processing ' , img , ' by >>> ' , ' (OMPvar,   CSMPSP,   gOMP,   ROMP,   StOMP,   NIHT,    HL0T,   IHT)\n'));
        
        alpha = 0.2;
        for ii = 1:1:sparse
            [OMPvar_ERROR(ii,1),OMPvar_SSIM(ii,1)]  = OMPvar_ano(Dic,Vector,ii,alpha,0);%2015
            [x_hat,CSMPSP_SSIM(ii,1),CSMPSP_ERROR(ii,1)] =   CSMPSP_ano(Vector,Dic,ii,{0 5});%2015
            [x_hat,gOMP_SSIM(ii,1),gOMP_ERROR(ii,1)]     =  gOMP_ano(Vector,Dic,ii,{0 5 4});%2012
            [x_hat,ROMP_SSIM(ii,1),ROMP_ERROR(ii,1)]    =   ROMP_ano(Vector,Dic,ii,{0 5});%2009
            [x_hat,StOMP_SSIM(ii,1),StOMP_ERROR(ii,1)]  =   StOMP_ano(Vector,Dic,ii,{0 5 3});%2012
            [x_hat,NIHT_SSIM(ii,1),NIHT_ERROR(ii,1)]    =   NIHT_ano(Vector,Dic,ii,{0 20});%2010
            [s, err_mse, iter_time,HL0T_SSIM(ii,1),HL0T_ERROR(ii,1)]  =  hard_l0_Mterm(Vector,Dic,size(Dic,2),ii);
            [xhat, r_mag,IHT_SSIM(ii,1),IHT_ERROR(ii,1)] = iht_ano(Vector, Dic, ii, 50, 0);
            check = 1;
        end
        RES = find(OMPvar_SSIM(:)>OMP_SSIM(:));
        RESU = find(OMPvar_ERROR(:)<OMP_ERROR(:));                                                
        fprintf(strcat('----------------------------,(%d%%, %d%%)'),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012
        RES = find(CSMPSP_SSIM(:)>OMP_SSIM(:));
        RESU = find(CSMPSP_ERROR(:)<OMP_ERROR(:));                                                        
        fprintf(strcat(',(%d%%, %d%%)'),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012   
        RES = find(gOMP_SSIM(:)>OMP_SSIM(:));
        RESU = find(gOMP_ERROR(:)<OMP_ERROR(:));                                                        
        
        fprintf(strcat(',(%d%%, %d%%)'),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012        
        RES = find(ROMP_SSIM(:)>OMP_SSIM(:));
        RESU = find(ROMP_ERROR(:)<OMP_ERROR(:));                                                        
        
        fprintf(strcat(',(%d%%, %d%%)'),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012        
        RES = find(StOMP_SSIM(:)>OMP_SSIM(:));
        RESU = find(StOMP_ERROR(:)<OMP_ERROR(:));                                                        
        
        fprintf(strcat(',(%d%%, %d%%)'),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012              
        RES = find(NIHT_SSIM(:)>OMP_SSIM(:));
        RESU = find(NIHT_ERROR(:)<OMP_ERROR(:));                                                        
        
        fprintf(strcat(',(%d%%, %d%%)'),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012        
        RES = find(HL0T_SSIM(:)>OMP_SSIM(:));
        RESU = find(NIHT_ERROR(:)<OMP_ERROR(:));                                                        
        
        fprintf(strcat(',(%d%%, %d%%) '),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012        
        RES = find(IHT_SSIM(:)>OMP_SSIM(:));
        RESU = find(IHT_ERROR(:)<OMP_ERROR(:));                                                        
        
        fprintf(strcat(',(%d%%, %d%%) \n'),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        fprintf(strcat('Processing ' , img , ' by >>> ' , ' (sOMP)'));
        [~, sOMP_ERROR, sOMP_SSIM, y_sOMP] = SSIM_selection_OMP_new(Dic, Vector,sparse);            
        RES = find(sOMP_SSIM(:)>OMP_SSIM(:));
        RESU = find(sOMP_ERROR(:)<OMP_ERROR(:));                                                               
        fprintf(strcat(',  (%d%%, %d%%)\n'),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        fprintf(strcat('Processing ' , img , ' by >>> ' , ' (GSM)'));
        [x_sol, sol, GSM_ERROR, GSM_SSIM] = GSM_ano(Dic,Vector,sparse,'profile', 'fast');%2021%'fast', 'normal' or 'thorough'
        RES = find(GSM_SSIM(:)>OMP_SSIM(:));
        RESU = find(GSM_ERROR(:)<OMP_ERROR(:));                                                                
        fprintf(strcat(',  (%d%%, %d%%)\n'),size(RES,1)*100/sparse,size(RESU,1)*100/sparse);%2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        SSIM_ALL = SSIM_ALL + [OMP_SSIM OMPvar_SSIM COSAMP_SSIM GSM_SSIM  sOMP_SSIM CSMPSP_SSIM gOMP_SSIM ROMP_SSIM StOMP_SSIM NIHT_SSIM SAMP_SSIM HL0T_SSIM];
        ERROR_ALL = ERROR_ALL  + [OMP_ERROR OMPvar_ERROR COSAMP_ERROR GSM_ERROR  sOMP_ERROR CSMPSP_ERROR gOMP_ERROR ROMP_ERROR StOMP_ERROR NIHT_ERROR SAMP_ERROR HL0T_ERROR];


        
        check = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(iOMP)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %[~, iomp_curERR, iomp_ssim, y_iOMP]                      = SSIM_iOMP_new(Dic, Vector,sparse);   
        %iomp_runtime = toc;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(MiOMP)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %tic;
        %[miomp_ssim, miomp_curERR, y_i_miomp]                      = getMIOMP(omp_ssim, iomp_ssim,OMP_curERR, iomp_curERR,y_i_omp,y_iOMP, sparse); 
        %comp_time = toc;  
        %miomp_time = omp_runtime + iomp_runtime + comp_time;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(nOMP)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %[~, coeff, iopt, OMP_curERR,...
        % nOMP_curERR,omp_ssim,nomp_ssim, y_i_omp, y_i_nomp,ssimoOPT_time]       = SSIM_OMP_nOMP_new(Dic, Vector,sparse);
        %nOMP_time = ssimoOPT_time + omp_runtime;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(sOMP)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %tic; 
        %[~, somp_curERR, somp_ssim, y_sOMP]                = SSIM_selection_OMP_new(Dic, Vector,sparse);             
        %somp_time = toc;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(LocBE)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %tic;
        %[~, LocBE_curERR,  LocBE_ssim, y_i_LocBE]                = LocBE_new(Dic, Vector,coeff,iopt,omp_ssim,y_i_omp);
        %LocBE_time = toc;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%(FBP)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %tic;
        %try [~, FBP_curERR, FBP_ssim, y_i_FBP]                       = FBP_new(Vector,Dic, 13,12,'K',0,sparse); 
        %catch me;
         %   FBP_ssim(1,1:sparse) = 0;
         %   FBP_curERR(1:sparse,1) = 0;
        %end
        %FBP_time = toc;
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %% temp
           %temp_comp_ssim = (temp_comp_ssim+[omp_ssim' somp_ssim']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %if size(FBP_ssim) < sparse
        %    FBP_ssim(1,1:sparse) = 0;
        %    FBP_curERR(1:sparse,1) = 0;
            
        %end
         y_i_omp_(:,:,i,j) =  y_i_omp;
         %y_i_iomp_(:,:,i,j) = y_iOMP;
         %y_i_nomp_(:,:,i,j) =  y_i_nomp;
        for kk =1:1:sparse
            [SS_omp_den(kk,:) NM_omp_den(kk,:) bm3d_runtime1] = getOptDenoised_OneStage(Vector, y_i_omp(:,kk),OMP_SSIM(kk),OMP_ERROR(kk));
            %[SS_miomp_den(kk,:) NM_miomp_den(kk,:) bm3d_runtime2] = getOptDenoised_OneStage(Vector, y_i_miomp(:,kk),miomp_ssim(kk),miomp_curERR(kk));
            %[SS_iomp_den(kk,1) NM_iomp_den(kk,1)] = getOptDenoised_OneStage(Vector, y_iOMP(:,kk),iomp_ssim(kk),iomp_curERR(kk));
            %[SS_nomp_den(kk,:) NM_nomp_den(kk,:) bm3d_runtime3] = getOptDenoised_OneStage(Vector, y_i_nomp(:,kk),nomp_ssim(kk),nOMP_curERR(kk));
            
            %[SS_somp_den(kk,:) NM_somp_den(kk,:)] = getOptDenoised_OneStage(Vector, y_sOMP(:,kk),somp_ssim(kk),somp_curERR(kk));
            %[SS_locbe_den(kk,:) NM_locbe_den(kk,:)] = getOptDenoised_OneStage(Vector, y_i_LocBE(:,kk),LocBE_ssim(kk),LocBE_curERR(kk));
  
        
        end
         %size_ = size(cell2mat((struct2cell(SS_omp_den)')));omp_runtime/sparse  iomp_runtime/sparse 
         %omp_t = omp_runtime/sparse;
         %omp_bm3d_t = omp_runtime/sparse+bm3d_runtime1(1);
         %miomp_t = miomp_time/sparse;
         %miomp_bm3d_t = miomp_time/sparse+bm3d_runtime1(1);
         %FBP_t = FBP_time/sparse;
         %LocBE_t = LocBE_time/sparse;
         
         ROS(:,:,k) =  SS_omp_den;
         RON(:,:,k) =  NM_omp_den;
         %RiS(:,k)   =  iomp_ssim;
         %RiN(:,k)   =  iomp_curERR; 
         %RFS(:,k)   =  FBP_ssim';
         %RFN(:,k)   =  FBP_curERR;  
         %RmS(:,:,k)   =  SS_miomp_den;
         %RmN(:,:,k)   =  NM_miomp_den;  
         
         %RnS(:,:,k) =  SS_nomp_den;
         %RnN(:,:,k) =  NM_nomp_den;
         %RsS(:,k)   =  somp_ssim'; 
         %RsN(:,k)   =  somp_curERR; 
         %RLS(:,k)   =  LocBE_ssim';
         %RLN(:,k)   =  LocBE_curERR;  
         
         
         ROS_m  =  ROS_m + ROS(:,:,k);
         RON_m  =  RON_m + RON(:,:,k);
         %RiS_m  =  RiS_m + RiS(:,k);
         %RiN_m  =  RiN_m + RiN(:,k);
         %RFS_m  =  RFS_m + RFS(:,k);
         %RFN_m  =  RFN_m + RFN(:,k);
         %RmS_m  =  RmS_m + RmS(:,:,k);
         %RmN_m  =  RmN_m + RmN(:,:,k);
         %RnS_m  =  RnS_m + RnS(:,:,k);
         %RnN_m  =  RnN_m + RnN(:,:,k);
         %RsS_m  =  RsS_m + RsS(:,k);
         %RsN_m  =  RsN_m + RsN(:,k);
         %RLS_m  =  RLS_m + RLS(:,k);
         %RLN_m  =  RLN_m + RLN(:,k);
         
         %COLLECTED_SSIM = COLLECTED_SSIM + [omp_ssim' iomp_ssim nomp_ssim somp_ssim'];
         %COLLECTED_NMSE = COLLECTED_NMSE + [OMP_curERR iomp_curERR nOMP_curERR somp_curERR];
         %COLLECTED_SSIM_avg = COLLECTED_SSIM./k;
         %COLLECTED_NMSE_avg = COLLECTED_NMSE./k;
         Check = 0;
             %RESULTS_FBP(k)= FBP(dic, Vector,coeff,iopt,omp_ssim);

             %vectorHat(i,j,:) = RESULTS_third(k).Xhat_OMP;
            %batchesHat(i,j,:,:)= vectorToBatch(vectorHat(i,j,:),batchSize);
            %fprintf('Processing batch %i : %i\n',index, int32(batchPerRow*batchPerCol));
          clc;
          fprintf('Processing %s >>> %.2f %%\n',img,index*100/(batchPerRow*batchPerCol));

          index = index + 1;  
          
          
          k = k + 1;
     end
end
         ROS_m  =  ROS_m./k;
         RON_m  =  RON_m./k;
         %RiS_m  =  RiS_m./k;
         %RiN_m  =  RiN_m./k;
         %RFS_m  =  RFS_m./k;
         %RFN_m  =  RFN_m./k;
         %RmS_m  =  RmS_m./k;
         %RmN_m  =  RmN_m./k;
         %RnS_m  =  RnS_m./k;
         %RnN_m  =  RnN_m./k;
         %RsS_m  =  RsS_m./k;
         %RsN_m  =  RsN_m./k;
         %RLS_m  =  RLS_m./k;
         %RLN_m  =  RLN_m./k;   
         check  = 0;
         RESULTS.ROS  =  ROS;
         RESULTS.RON  =  RON;
         %RESULTS.RiS  =  RiS;
         %RESULTS.RiN  =  RiN; 
         %RESULTS.RFS  =  RFS;
         %RESULTS.RFN  =  RFN;  
         %RESULTS.RmS  =  RmS;
         %RESULTS.RmN  =  RmN;
         %RESULTS.RnS  =  RnS;
         %RESULTS.RnN  =  RnN;
         %RESULTS.RsS  =  RsS; 
         %RESULTS.RsN  =  RsN; 
         %RESULTS.RLS  =  RLS;
         %RESULTS.RLN  =  RLN;  

         RESULTS.ROS_m  =  ROS_m;
         RESULTS.RON_m  =  RON_m;
         %RESULTS.RiS_m  =  RiS_m;
         %RESULTS.RiN_m  =  RiN_m; 
         %RESULTS.RFS_m  =  RFS_m;
         %RESULTS.RFN_m  =  RFN_m;  
         %RESULTS.RmS_m  =  RmS_m;
         %RESULTS.RmN_m  =  RmN_m;
         %RESULTS.RnS_m  =  RnS_m;
         %RESULTS.RnN_m  =  RnN_m;
         %RESULTS.RsS_m  =  RsS_m; 
         %RESULTS.RsN_m  =  RsN_m; 
         %RESULTS.RLS_m  =  RLS_m;
         %RESULTS.RLN_m  =  RLN_m; 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%% we do the for_back processing for the main algs%%%%%%%%
         %%%%%%%%% OMP, iOMP, FBP, MiOMP, nOMP, LocBE, sOMP%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %[SS_den, NM_den, imgs] = buildImages(y_i_omp_,y_i_iomp_,y_i_nomp_,floor(batchPerRow),floor(batchPerCol),img_pixels,img_double);

         %RESULTS.FOR_RESULTS = forProcessing(RESULTS, sparse,img);
         %RESULTS.BAC_RESULTS = bacProcessing(RESULTS, sparse,img);
         
         check = 0;

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










