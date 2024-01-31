%good for nomp
%(1039,1033,1031,1027,1017,1016,1015,1014,1013,1011,1010,1007,1006,1001 semi)
%(2039 semi,2033,2027,2001 semi
%(20,30, 40, 50,120,130,140,160,210,250,270,290,310 semi,360,400,410,460,
%semi,520, 550, 560,570,600 semi, 610,620,640 semi, 670,740 semi,750,760
%semi, 770, 780 semi, 800,810, 820,830,840, 850,880, 940,950,900 semi,930 vg,970,980, 990 vg, 1000 semi
%bad for all (1030,180,890, 
function  runCollect_BESIR()
clc;
count = 1;
for i=1:1:1000%10:10:1000
    t1 = clock;
    
    str1 = 'AAA__file_(';%'BBB__file_('
    str2 = int2str(i);
    str3 = ').mat';
    matFile = strcat(str1,str2,str3);
    fprintf('Processing file ( %s ) >>> ',matFile);

    %matFile = 'BBB__file_(1006).mat';%I stopped at 510
    %RES{count}.res = collectResults_BESIR(matFile);
    %RES{count}.res = collectResults_BESIR_newAlgs(matFile);
    RES(i).res = collectResults_MiOMP(matFile);
    %res = collectResults_MiOMP(matFile);
    %save_(res,i);
    t2 = clock;
    fprintf(' >>> %f minutes \n',etime(t2,t1)/60);
    %fprintf(' >>> No. Files %i\n',count);
    %RES1 = collectResults_BESIR_inc_SOMP(matFile);
    count = count + 1;
    x = 0;
end
save('AAA_MiOMP_1000_128_128','RES');
%save('AAA_MiOMP_40_128_128','RES');
%save('MiOMP_1000_128_128','RES');
%save('BIG_DATASET_RESULTS','RES');
xx = 0;

end
function save_(RES,i)
    save(strcat('AAA_MiOMP_1000_128_128_', int2str(i)),'RES');

end