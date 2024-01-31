function  runAVGs_MultipleFiles()
Lear = 'AAA_sp_file_(';
Stru = 'BBB_sp_file_(';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROS = []; RiS = []; RFS = []; RmS = []; RnS = []; RLS = []; RsS = [];
RON = []; RiN = []; RFN = []; RmN = []; RnN = []; RLN = []; RsN = [];
h = 1;
%for file = 220001:1:220048%450004:1:450004; 650001:1:650004
%for file = 190001:1:190024   
%%%%for file = 210001:1:210048% except 210009   
%for file = 310001:1:310100%    
%for file = 330001:1:330100%  
%for file = 360001:1:360100% 
%for file = 410001:1:410064%    
%%%%%for file = 420001:1:420015%  except 420015  upto 420064  
%for file = 430001:1:430064% 
%for file = 460001:1:460064% 
%for file = 510001:1:510040%
%%%%%%%%%%for file = 520001:1:520040%
For_RES = [];
RO_ALL = [];
RI_ALL = [];
RF_ALL = [];
RM_ALL = [];
RN_ALL = [];
RS_ALL = [];
RL_ALL = [];

b = 1:1:12;
k = b';
%for file = 560001:1:560040%     
%for file = 650001:1:650001%     
for file = 250003:1:250003 %190001:1:190024    
    d = int2str(file);
    %fi = strcat(Lear,d,')');
    fi = strcat(Stru,d,')');
    x = load(fi);%RES_MiOMP_AAA_DS_1000_64; %RES_BBB_DS_128;%RES_BBB_DS_256;
    %x.RESULTS.ROS; 
    [A, B, C] = size(x.RESULTS.ROS);
    fin = ((h-1)*C+1);
    lin = h*C;
    disp(fi);
    ROS(:,:,fin:lin) = x.RESULTS.ROS; 
    RiS(:,fin:lin) = x.RESULTS.RiS; 
    RFS(:,fin:lin) = x.RESULTS.RFS; 
    RmS(:,:,fin:lin) = x.RESULTS.RmS; 
    RnS(:,:,fin:lin) = x.RESULTS.RnS; 
    RLS(:,fin:lin) = x.RESULTS.RLS; 
    RsS(:,fin:lin) = x.RESULTS.RsS; 
    RON(:,:,fin:lin) = x.RESULTS.RON; 
    RiN(:,fin:lin) = x.RESULTS.RiN; 
    RFN(:,fin:lin) = x.RESULTS.RFN; 
    RmN(:,:,fin:lin) = x.RESULTS.RmN; 
    RnN(:,:,fin:lin) = x.RESULTS.RnN; 
    RLN(:,fin:lin) = x.RESULTS.RLN; 
    RsN(:,fin:lin) = x.RESULTS.RsN; 
    for ii = fin:1:lin
        RO_ALL = [RO_ALL; k ROS(:,:,ii) RON(:,:,ii)];
        RI_ALL = [RI_ALL; k RiS(:,ii) RiN(:,ii)];
        RF_ALL = [RF_ALL; k RFS(:,ii) RFN(:,ii)];
        RM_ALL = [RM_ALL; k RmS(:,:,ii) RmN(:,:,ii)];
        RN_ALL = [RN_ALL; k RnS(:,:,ii) RnN(:,:,ii)];
        RS_ALL = [RS_ALL; k RsS(:,ii) RsN(:,ii)];
        RL_ALL = [RL_ALL; k RLS(:,ii) RLN(:,ii)];
    end
    %%%%%
    For_RES = [For_RES; x.RESULTS.FOR_RESULTS.For_RESULTS_SS_all];
    h = h + 1;
    %ROS_ALL = [ROS_ALL; x.RESULTS.ROS];
end
R_ALL = [RO_ALL RI_ALL RF_ALL RM_ALL RN_ALL RS_ALL RL_ALL];
R_ALL_mean = [];
for j=1:1:12
    m = j:12:size(R_ALL,1);
        R_ALL_mean = [R_ALL_mean; nanmean(R_ALL(m,:))];
        cc=0;
end
    
ros_rho = [];
ris_rho = [];
rfs_rho = [];
rms_rho = [];
rns_rho = [];
rls_rho = [];
rss_rho = [];
ros_m = [];
ris_m = [];
rfs_m = [];
rms_m = [];
rns_m = [];
rls_m = [];
rss_m = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ron_rho = [];
rin_rho = [];
rfn_rho = [];
rmn_rho = [];
rnn_rho = [];
rln_rho = [];
rsn_rho = [];
ron_m = [];
rin_m = [];
rfn_m = [];
rmn_m = [];
rnn_m = [];
rln_m = [];
rsn_m = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A, B, C] = size(ROS);
for i = 1:1:A
              clc;
          fprintf('Processing %.2f, %.2f >>>',i,A);
    k = 1;
    for j = 1:1:C
        ros(j,:)  = ROS(i,:,j);
        ris(j,:) = RiS(i,j);
        rfs(j,:) = RFS(i,j);
        rms(j,:) = RmS(i,:,j);
        rns(j,:) = RnS(i,:,j);
        rls(j,:) = RLS(i,j);
        rss(j,:) = RsS(i,j);
        ron(j,:)  = RON(i,:,j);
        rin(j,:) = RiN(i,j);
        rfn(j,:) = RFN(i,j);
        rmn(j,:) = RmN(i,:,j);
        rnn(j,:) = RnN(i,:,j);
        rln(j,:) = RLN(i,j);
        rsn(j,:) = RsN(i,j);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    rss_m = [rss_m; nanmean(rss)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ron_rho = [ron_rho; my_ttest(ron(:,1), ron(:,2)) my_ttest(ron(:,1), ron(:,3)) my_ttest(ron(:,1), ron(:,4)) my_ttest(ron(:,1), ron(:,5)) my_ttest(ron(:,1), ron(:,6))];
    rin_rho = [rin_rho; my_ttest(ron(:,1), rin(:,1))];
    rfn_rho = [rfn_rho; my_ttest(ron(:,1), rfn(:,1))];
    rmn_rho = [rmn_rho; my_ttest(ron(:,1), rmn(:,1)) my_ttest(ron(:,1), rmn(:,2)) my_ttest(ron(:,1), rmn(:,3)) my_ttest(ron(:,1), rmn(:,4)) my_ttest(ron(:,1), rmn(:,5)) my_ttest(ron(:,1), rmn(:,6))];
    rnn_rho = [rnn_rho; my_ttest(ron(:,1), rnn(:,1)) my_ttest(ron(:,1), rnn(:,2)) my_ttest(ron(:,1), rnn(:,3)) my_ttest(ron(:,1), rnn(:,4)) my_ttest(ron(:,1), rnn(:,5)) my_ttest(ron(:,1), rnn(:,6))];
    rln_rho = [rln_rho; my_ttest(ron(:,1), rln(:,1))];
    rsn_rho = [rsn_rho; my_ttest(ron(:,1), rsn(:,1))];
    ron_m = [ron_m; nanmean(ron)];
    rin_m = [rin_m; nanmean(rin)];
    rfn_m = [rfn_m; nanmean(rfn)];
    rmn_m = [rmn_m; nanmean(rmn)];
    rnn_m = [rnn_m; nanmean(rnn)];
    rln_m = [rln_m; nanmean(rln)];
    rsn_m = [rsn_m; nanmean(rsn)];    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cc = 0;
    ALL_nOMP_paper = [ros_m(:,1) ris_m(:,1) rns_m(:,1) ros_m(:,2) rns_m(:,2)];
    ALL_MiOMP_paper = [ros_m(:,1) ris_m(:,1) rms_m(:,1) ros_m(:,2) rms_m(:,2)];
    ALL_LocBE_paper = [ros_m(:,1) ros_m(:,2) ros_m(:,3) ros_m(:,4) ros_m(:,5) ros_m(:,6) rfs_m(:,1) rls_m(:,1)];
    ALL_LocBE_pap3 = [ros_m(:,1) ros_m(:,2) ros_m(:,3) ros_m(:,4) ros_m(:,5) rms_m(:,1) rms_m(:,2) rfs_m(:,1) rls_m(:,1)];
    ALL_sOMP_paper = [ros_m(:,1) ros_m(:,2) rss_m(:,1)];
    cust_ss_pap3 = [ros_m(:,1) ros_m(:,2)  rms_m(:,1) rms_m(:,2) rfs_m(:,1) rls_m(:,1)];

    
    ALL_rho_nOMP = [ris_rho(:,1) rns_rho(:,1) ros_rho(:,1) rns_rho(:,2)];
    ALL_rho_MiOMP = [ris_rho(:,1) rms_rho(:,1) ros_rho(:,1) rms_rho(:,2)];
    ALL_rho_LocBE = [ros_rho(:,1) ros_rho(:,2) ros_rho(:,3) ros_rho(:,4) ros_rho(:,5) rfs_rho(:,1) rls_rho(:,1)];
    ALL_rho_LocBE_pap3 = [ros_rho(:,1) ros_rho(:,2) ros_rho(:,3) ros_rho(:,4) ros_rho(:,5) rms_rho(:,1) rms_rho(:,2) rfs_rho(:,1) rls_rho(:,1)];
    ALL_rho_sOMP = [rfs_rho(:,1) rss_rho(:,1)];
    cust_rho_pap3 = [ros_rho(:,1) rms_rho(:,1) rms_rho(:,2) rfs_rho(:,1) rls_rho(:,1)];
   %BM3D TVL1 SBATV BILATERAL NLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ALL_nOMP_paper_n = [ron_m(:,1) rin_m(:,1) rnn_m(:,1) ron_m(:,2) rnn_m(:,2)];
    ALL_MiOMP_paper_n = [ron_m(:,1) rin_m(:,1) rmn_m(:,1) ron_m(:,2) rmn_m(:,2)];
    ALL_LocBE_paper_n = [ron_m(:,1) ron_m(:,2) rfn_m(:,1) rln_m(:,1)];
    ALL_sOMP_paper_n = [ron_m(:,1) ron_m(:,2) rsn_m(:,1)];

    ALL_rho_nOMP_n = [rin_rho(:,1) rnn_rho(:,1) ron_rho(:,1) rnn_rho(:,2)];
    ALL_rho_MiOMP_n = [rin_rho(:,1) rmn_rho(:,1) ron_rho(:,1) rmn_rho(:,2)];
    ALL_rho_LocBE_n = [rfn_rho(:,1) rln_rho(:,1)];
    ALL_rho_sOMP_n = [rsn_rho(:,1)];

    %ros_all = [ros_all;ros] ;
    
 
end
y=0;
%x.x;
GI_sel = 1:1:48;%GI_7;%1:40;%GI_5;%1:40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end