clc;

file = load('AAA_MiOMP_1000_32_32');
res = file.RES;
size = length(res);
count = 1;
temp = [];
for i=1:1:size
    RES_ = res(i).res.SSIM_Fixed_sp_res;
    ratio = (RES_(:,4)-RES_(:,2))*100./RES_(:,2);
    xx = 0;
    temp = [temp ratio];
end
cc = 0;
