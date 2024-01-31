clc;
%load('BBB_sp_file_(220046).mat')
%load('BBB_sp_file_(220047).mat')
%load('BBB_sp_file_(520025).mat')
    %fi = strcat(Lear,d,')');
    
init = 250003;
 fin = 250003;
 count = 0;
 count_locbe = 0;
 count_bm3d = 0;
 count_miomp = 0;
 count_miomp_bm3d = 0;
 count_fbp = 0;
 prob_all = [];
for iii =  init:1:fin
    d = int2str(iii);
    fi = strcat('BBB_sp_file_(',d,').mat');
    file = load(fi);
    res = file.RESULTS;
    ROS = res.ROS;
    RmS = res.RmS;
    RON = res.RON;
    RFS = res.RFS;
    RLS = res.RLS;
    RLN = res.RLN;
    size1 = size(ROS,1);
    size3 = size(ROS,3);
    
    temp = [];
    prob = [];
    for i=1:1:size1
        for j=1:1:size3
            if RLS(i,j)>ROS(i,1,j)
                count = count + 1;
                count_locbe = count_locbe + 1;
            end
            if ROS(i,2,j)>ROS(i,1,j)
                %count = count + 1;
                count_bm3d = count_bm3d + 1;
            end
            if RmS(i,1,j)>ROS(i,1,j)
                %count = count + 1;
                count_miomp = count_miomp + 1;
            end   
            if RmS(i,2,j)>ROS(i,1,j)
                %count = count + 1;
                count_miomp_bm3d = count_miomp_bm3d + 1;
            end
            if RFS(i,j)>ROS(i,1,j)
                %count = count + 1;
                count_fbp = count_fbp + 1;
            end
        end
     prob = [prob;count/size3];
        count = 0;
    end
    prob_omp_bm3d = count_bm3d/(size3*12);
    prob_miomp = count_miomp/(size3*12);
    prob_miomp_bm3d = count_miomp_bm3d/(size3*12);
    prob_fbp = count_fbp/(size3*12);
    prob_locbe = count_locbe/(size3*12);
    prob_all = [prob_all;prob_omp_bm3d prob_miomp prob_miomp_bm3d prob_fbp prob_locbe]; 
    count = 0;
    count_locbe = 0;
    count_bm3d = 0;
    count_miomp = 0;
    count_miomp_bm3d = 0;
    count_fbp = 0; 
end
prob_omp_bm3d = count_bm3d/(size3*(fin-init+1)*12);
prob_miomp = count_miomp/(size3*(fin-init+1)*12);
prob_miomp_bm3d = count_miomp_bm3d/(size3*(fin-init+1)*12);
prob_fbp = count_fbp/(size3*(fin-init+1)*12);
prob_locbe = count_locbe/(size3*(fin-init+1)*12);
prob_all = [prob_omp_bm3d prob_miomp prob_miomp_bm3d prob_fbp prob_locbe];
cc = 0;
