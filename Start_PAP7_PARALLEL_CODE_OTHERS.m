aa = [320001];
for i = 320001:1:322000
    try
        zz = load(strcat('PAP7_OTHERS_file_(',int2str(i),').mat'));
    catch
        aa = [aa,i];
    end
end




addpath(genpath ('/home/matlabuser2'));
%aa = [351:1:359, 434:1:446, 459:1:462, 567:1:576, 735:1:757, 927:1:929, 951,1011:1:1021,1047,1113:1:1118,1149:1:1154,1185:1:1196];
%aaa = 320000 + aa;
%parfor file_index = 320001:322000
for file_index = 1:size(aa,2)

   % Start_PAP7_OTHERS(file_index);
    Start_PAP7_OTHERS(aa(1,file_index));
end