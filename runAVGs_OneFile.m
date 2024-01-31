function  runAVGs_OneFile()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = load('AAA_sp_file_(650001)');%RES_MiOMP_AAA_DS_1000_64; %RES_BBB_DS_128;%RES_BBB_DS_256;
ROS = x.RESULTS.ROS;    ros_rho = [];
[A, B, C] = size(ROS);
RiS = x.RESULTS.RiS; ris_rho = [];
RFS = x.RESULTS.RFS; rfs_rho = [];
RmS = x.RESULTS.RmS; rms_rho = [];
RnS = x.RESULTS.RnS; rns_rho = [];
RLS = x.RESULTS.RLS; rls_rho = [];
RsS = x.RESULTS.RsS; rss_rho = [];
ros_m = [];
ris_m = [];
rfs_m = [];
rms_m = [];
rns_m = [];
rls_m = [];
rss_m = [];
for i = 1:1:A
    k = 1;
    for j = 1:1:C
        ros(j,:)  = ROS(i,:,j);
        ris(j,:) = RiS(i,j);
        rfs(j,:) = RFS(i,j);
        rms(j,:) = RmS(i,:,j);
        rns(j,:) = RnS(i,:,j);
        rls(j,:) = RLS(i,j);
        rss(j,:) = RsS(i,j);
    end
    ros_rho = [ros_rho; my_ttest(ros(:,1), ros(:,2)) my_ttest(ros(:,1), ros(:,3)) my_ttest(ros(:,1), ros(:,4)) my_ttest(ros(:,1), ros(:,5)) my_ttest(ros(:,1), ros(:,6))];
    ris_rho = [ris_rho; my_ttest(ros(:,1), ris(:,1))];
    rfs_rho = [rfs_rho; my_ttest(ros(:,1), rfs(:,1))];
    rms_rho = [rms_rho; my_ttest(ros(:,1), rms(:,1)) my_ttest(ros(:,1), rms(:,2)) my_ttest(ros(:,1), rms(:,3)) my_ttest(ros(:,1), rms(:,4)) my_ttest(ros(:,1), rms(:,5)) my_ttest(ros(:,1), rms(:,6))];
    rns_rho = [rns_rho; my_ttest(ros(:,1), rns(:,1)) my_ttest(ros(:,1), rns(:,2)) my_ttest(ros(:,1), rns(:,3)) my_ttest(ros(:,1), rns(:,4)) my_ttest(ros(:,1), rns(:,5)) my_ttest(ros(:,1), rns(:,6))];
    rls_rho = [rls_rho; my_ttest(ros(:,1), rls(:,1))];
    rss_rho = [rss_rho; my_ttest(ros(:,1), rss(:,1))];
    ros_m = [ros_m; nanmean(ros)];
    ris_m = [ris_m; nanmean(ris)];
    rfs_m = [rfs_m; nanmean(rfs)];
    rms_m = [rms_m; nanmean(rms)];
    rns_m = [rns_m; nanmean(rns)];
    rls_m = [rls_m; nanmean(rls)];
    rss_m = [rss_m; nanmean(rss)];    cc = 0;
    %ros_all = [ros_all;ros] ;
    
 
end
y=0;
%x.x;
GI_sel = 1:1:48;%GI_7;%1:40;%GI_5;%1:40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end