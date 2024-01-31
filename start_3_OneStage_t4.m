%imDim,imName,imExt,dicType,ssimDIF, enh
sparse = 12;
imDim = 32; %512,256,128
imExt ='.jpg'; %jpg, gif
DATASET(1,1) = 190001 ;DATASET(1,2) = 190024;
DATASET(2,1) = 250001 ;DATASET(2,2) = 250048;
DATASET(3,1) = 220001 ;DATASET(3,2) = 220048;
DATASET(4,1) = 210001 ;DATASET(4,2) = 210048;
DATASET(5,1) = 310001 ;DATASET(5,2) = 310100;
DATASET(6,1) = 360001 ;DATASET(6,2) = 360100;%start
DATASET(7,1) = 330001 ;DATASET(7,2) = 330100;%start
DATASET(8,1) = 420001 ;DATASET(8,2) = 420016;
DATASET(9,1) = 5200011 ;DATASET(9,2) = 520040;% wait
DATASET(10,1) = 510001 ;DATASET(10,2) = 510040;
DATASET(11,1) = 560001 ;DATASET(11,2) = 560040;
DATASET(12,1) = 410001 ;DATASET(12,2) = 410064;
DATASET(13,1) = 460001 ;DATASET(13,2) = 460064;
DATASET(14,1) = 430001 ;DATASET(14,2) = 430064;
c = [650001, 450004];

for j = 14:1:14
for i = DATASET(j,1):1:DATASET(j,1)%120:1:1000%753:5:998
    t1 = clock;
    str1 = '(';%1,15,16,35,47,48 v
    str2 = int2str(i);
    str3 = ')';
    imName = strcat(str1,str2,str3);
    % .gif, .tif
    dicType = 'learned'; %learned, structured
    ssimDIF = 0;
    enh = 0; %0, 1 to use the ssim enhancement idea
    
    stopCRITERIA = 1;% 0 - (L2 error), 1 - (SPARSE),2 - (ssim), 
    quantize_bits = 8;
    %fprintf('Processing file ( %s ) >>> ',imName);
    doMyFunction_BE_OneStage(imDim,imName,imExt,dicType,ssimDIF, enh, sparse,quantize_bits);
    t2 = clock;
    %fprintf('... >>> completed in  %f minutes \n',etime(t2,t1)/60);

    cc = 0;
end
end