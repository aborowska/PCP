clear all
close all

model = 'ar1';
sigma1 = 1;
sigma2 = 2;
H = 100;
v_new = '(R2015a)';
II = 100;

T = 100;
DM_100 = DM_test(model, T, sigma1, sigma2, H, II);

T = 1000;
DM_1000 = DM_test(model, T, sigma1, sigma2, H, II);

T = 10000;
DM_10000 = DM_test(model, T, sigma1, sigma2, H, II);


