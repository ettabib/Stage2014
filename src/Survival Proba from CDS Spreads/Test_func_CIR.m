clear all
clc

param = [0.143; 0.534];

s = 0:0.01:10; 


phi_1 = phi_CIR(s, param);
phi_2 = phi_CIR_2(s, param);
 
plot(s, phi_1);
hold on;
plot(s, phi_2, 'r');

xi_1 = xi_CIR(s, param);
xi_2 = xi_CIR_2(s, param);

figure;
plot(s, xi_1);
hold on;
plot(s, xi_2, 'r');
