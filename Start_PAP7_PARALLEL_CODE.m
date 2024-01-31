
%save('file_index.mat',"file_index");
%file_ = load('file_index.mat');
%file_index = file_.file_index;
%addpath(genpath ('/home/matlabuser2'));
%file_index = [220001:1:220048, 321001:1:322000];
%file_index_ = file_index(1:24);

%ss = size(file_index,2);
%numGPUs = gpuDeviceCount("available");
%parpool(numGPUs);
%garray = gpuArray(file_index_);
parfor s = 310801:311000%310101:310800%210001:210047%310801:311000 %310101:310600%210001:210047%310101:310300%321001:321200%321001:321500%220001:220048%  321001:321500


    %Start_PAP7_CODE(file_index(index));
    %Start_PAP7_OMP_EOMP_CODE_extra(file_index(index));
    %Start_PAP7_OOMP_EOOMP_CODE_extra(file_index(index));
    %Start_PAP7_FBP_EFBP_CODE_extra(file_index(index));
    %Start_PAP7_OMP_EOMP_CODE_extra_var1(s);
    %Start_PAP7_OMP_EOMP_CODE_extra_var2(s);
    %Start_PAP7_OMP_EOMP_CODE_extra_var3(s);
    %Start_PAP7_SISM_v_4(s);
     %Start_PAP7_SISM_v_4_All_SelfTaught(s);


     %Start_PAP7_SISM_v_4_EXP2_get_Dics_128(s)
     %Start_PAP7_SISM_v_4_EXP2_get_Dics(s);
     %Start_PAP7_SISM_v_4_EXPs_3_DS2(s)



    %Start_PAP7_SISM_v_4_EXPs_1_2_FINAL_FINAL(s)
     Start_PAP7_SISM_v_4_EXPs_4_DS2_FINAL_FINAL(s)
end