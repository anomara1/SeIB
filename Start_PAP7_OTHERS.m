%imDim,imName,imExt,dicType,ssimDIF, enh
function Start_PAP7_OTHERS(file_index)
clc;
batchSize = 256;%16,8
sparse = 90;
imDim = 32; %512,256,128
imExt ='.jpg'; %jpg, gif,jpeg


warning off;
for j = 6:1:6
    for i = file_index:file_index%DATASET(j,1):1:DATASET(j,2)%120:1:1000%753:5:998
        %t1 = clock;
        str1 = '(';%1,15,16,35,47,48 v
        str2 = int2str(i);
        str3 = ')';
        imName = strcat(str1,str2,str3);
        jpg_qual = [0:1:100];
        ssimDIF = 0;
        enh = 0;
        stopCRITERIA = 1;
        quantize_bits = 8;
        AAA_doMyFunction_BE_OneStage_PAP5__Others(batchSize,imDim,imName,imExt,ssimDIF, enh, sparse,quantize_bits);
        %AAA_doMyFunction_BE_OneStage_FBP_PAP5(batchSize,imDim,imName,imExt,ssimDIF, enh, sparse,quantize_bits);
    end
end
end
