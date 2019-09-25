close all;
clear;
clc;


T = readtable('1Blade_3D_VV_comparison_different_sides_PEC_mat_er_5_10GHz.dat');

Phi = T.PlaneWavePhi_deg_;

%T6 = readtable('dipole_FF_300e6_numseg31.dat');
figure;
plot(Phi, T.x1Bladde3D_all_er5_dBsm_, 'LineWidth', 2);
hold on;
plot(Phi, T.x1Bladde3D_Plastic_V2_er7_oneside263PEC_dBsm_, '-.', 'LineWidth', 2);
hold on;
plot(Phi, T.x1Bladde3D_Plastic_V2_er7_oneside261PEC_dBsm_, '*', 'LineWidth', 1);
hold on;
plot(Phi, T.x1Bladde3D_Plastic_V2_er7_oneside251_259_260PEC_dBsm_, '-.', 'LineWidth', 2);
hold on;
plot(Phi, T.x1Bladde3D_Plastic_V2_er7_oneside250PEC_dBsm_, '*', 'LineWidth', 1);
hold on;
plot(Phi, T.x1Bladde3D_Plastic_V2_er7_oneside256PEC_dBsm_, '-.', 'LineWidth', 2);
hold on;
grid on;

xlabel('Plane wave \phi (deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('RCS [dbm^2] at 3 GHz for 1 blade', 'FontSize', 12, 'FontWeight', 'bold');
title('RCS of 1 Blade 3D with different sides PEC', 'FontSize', 12, 'FontWeight', 'bold');

legend({'All \epsilon_r = 5', 'AB PEC', 'BC PEC', 'CD PEC', 'DE PEC', 'EA PEC'}, 'FontSize', 12, 'FontWeight', 'bold');

print('1Blades_RCS_comparison_different_side_PEC_VV_10GHz', '-depsc');
% %% Imp
% 
% close all;
% clear;
% clc;
% 
% 
% T = readtable('A3_A1_imp_5.dat');
% T1 = readtable('A3_A1_imp_61.dat');
% 
% %T6 = readtable('dipole_FF_300e6_numseg31.dat');
% figure;
% plot(T.Var1, T.Var2, 'LineWidth', 2);
% hold on;
% plot(T.Var1, T.Var3, 'LineWidth', 2);
% hold on;
% plot(T1.Var1, T1.Var2, '*', 'LineWidth', 2);
% hold on;
% plot(T1.Var1, T1.Var3, '-.', 'LineWidth', 2);
% hold on;
% grid on;
% 
% xlabel('Frequency (MHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('|Z_{in}| \Omega', 'FontSize', 12, 'FontWeight', 'bold');
% title('Impedance vs frequency', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Numseg = 5 With A3', 'Numseg = 5 with A1', 'Numseg = 61 with A3', 'Numseg = 61 with A1'}, 'FontSize', 12, 'FontWeight', 'bold');
% 
