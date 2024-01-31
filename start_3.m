%imDim,imName,imExt,dicType,ssimDIF, enh
imDim = 256; %512,256,128
imExt ='.gif'; %jpg, gif
for i = 2048:1:2048%1280050:1:1280100%120:1:1000%753:5:998
    t1 = clock;
    str1 = '(';%1,15,16,35,47,48
    str2 = int2str(i);
    str3 = ')';
    imName = strcat(str1,str2,str3);
    % .gif, .tif
    dicType = 'structured'; %learned, structured
    ssimDIF = 0;
    enh = 0; %0, 1 to use the ssim enhancement idea
    
    stopCRITERIA = 1;% 0 - (L2 error), 1 - (SPARSE),2 - (ssim), 
    quantize_bits = 8;
    %fprintf('Processing file ( %s ) >>> ',imName);
    doMyFunction_BE(imDim,imName,imExt,dicType,ssimDIF, enh, stopCRITERIA,quantize_bits);
    t2 = clock;
    %fprintf('... >>> completed in  %f minutes \n',etime(t2,t1)/60);

    cc = 0;
end
