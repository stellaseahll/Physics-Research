clear;clc;clf;

[t meanX meanP meanN meanE] = AlickiSpringNoise(50,10000,1,0,1,1,0,1);
save('Spring_T50_Nt100000_nh1_nc0_w1_r1_k1.mat');
[t meanX meanP meanN meanE] = AlickiSpringNoise(50,10000,1,0,1,0.1,0,1);
save('Spring_T50_Nt100000_nh1_nc0_w1_r0.1_k1.mat');
[t meanX meanP meanN meanE] = AlickiSpringNoise(50,10000,1,0,1,10,0,1);
save('Spring_T50_Nt100000_nh1_nc0_w1_r10_k1.mat');
[t meanX meanP meanN meanE] = AlickiSpringNoise(50,10000,1,0,1,1,0,0.1);
save('Spring_T50_Nt100000_nh1_nc0_w1_r1_k0.1.mat');
[t meanX meanP meanN meanE] = AlickiSpringNoise(50,10000,1,0,1,1,0,10);
save('Spring_T50_Nt100000_nh1_nc0_w1_r1_k10.mat');