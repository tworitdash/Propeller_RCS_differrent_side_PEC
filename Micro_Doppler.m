%% with different frequencies


close all;
clear;
clc;

f = [3.315e9, 10e9]; % place the frequencies used in FEKO. Size should be equal to size(ephi, 2) [Line number 57]
c = 3e8;
lamb = c./f;

% 4 GHz 
%filename = '../08202019_feko_one_prop_1_blade_Free_Space.out';
filename = '../08202019_feko_one_prop_2_blades_PEC.out';
%filename = '../08202019_feko_one_prop_1_blade_er2.out';
%filename = '../08202019_feko_one_prop_1_blade_er11.out';
%filename = '../08202019_feko_one_prop_1_blade_er81.out';
%filename = '../08202019_feko_one_prop_1_blade_er1_medium_dieled.out';
%filename = '../08202019_feko_one_prop_2_blades_er1_carbide.out';
%filename = '../08202019_Feko_two_blades_PEC.out';

efeReader_Feko(filename);

%fn = fullfile('../', '08202019_feko_one_prop_1_blade_Free_Space_data.txt');
fn = fullfile('../', '08202019_feko_one_prop_2_blades_PEC_data.txt');
%fn = fullfile('../', '08202019_feko_one_prop_1_blade_er2_data.txt');
%fn = fullfile('../', '08202019_feko_one_prop_1_blade_er11_data.txt');
%fn =  fullfile('../', '08202019_feko_one_prop_1_blade_er81_data.txt');
%fn = fullfile('../', '08202019_feko_one_prop_1_blade_er1_medium_dieled_data.txt');
%fn = fullfile('../', '08202019_feko_one_prop_2_blades_er1_carbide_data.txt');
%fn = fullfile('../', '08202019_Feko_one_blade_PEC_data.txt');
result_file = fn;
data = importdata(result_file);
data = regexp(strtrim(data), '\s+', 'split');
%data = reshape(data, 361, length(data)/361);


for data_i = 1:size(data,1)
    phi_angle_2(data_i) = str2num(data{data_i,1}{2}); % azimuth angle [degree]
    %E_theta_2(data_i) = (str2num(data{data_i,1}{3}).*str2num(data{data_i,1}{4})); % E_theta
    E_phi_2(data_i) = (str2num(data{data_i,1}{5}).* exp(1j * str2num(data{data_i,1}{6}) .* pi/180)); % E_phi
    
    rcs_sqm_2(data_i) = str2num(data{data_i,1}{7}); % RCS [m^2]
end



phi_angle_ = reshape(phi_angle_2, 10001, length(data)/10001);
%ephi_db_2 = 10*log10(E_phi_2);
ephi = reshape(E_phi_2, 10001, length(E_phi_2)/10001);

% phi_angle_ = reshape(phi_angle_2, 1441, length(data)/1441);
 rcs_db_2 = 10*log10(rcs_sqm_2);
 rcs = reshape(rcs_db_2, 10001, length(rcs_db_2)/10001);
% rcs_db_ = reshape(rcs_db_2, 1441, length(rcs_db_2)/1441);
%% 

close all;
for i = 1:size(ephi, 2)

rcs_i = rcs(:, i);
    
figure(i);
hold on;
plot(phi_angle_(:, 1), rcs_i, 'LineWidth', 2, 'color', [0.6350, 0.0780, 0.1840]);
grid on;
hold on;

xlabel('Plane wave \phi [deg]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('RCS [dBm^2]', 'FontSize', 12, 'FontWeight', 'bold');
title('RCS of the propeller', 'FontSize', 12, 'FontWeight', 'bold');
%legend({'Carbide', 'PEC', 'Free Space'}, 'FontSize', 12, 'FontWeight', 'bold');



ephi_i = ephi(:, i);
%Micro Doppler

wlen = 1024;                        % window length (recomended to be power of 2)
hop = wlen/8*7;                       % hop size (recomended to be power of 2)
nfft = wlen;                        % number of fft points (recomended to be power of 2)

% Calculating Fs from PRI
Omega = 6000; % rpm
Omega_s = Omega/60; % rps



L = 21.7e-2; %length of the blades
V_tip = 2 * pi * L * Omega_s;

f_dmax = (2 * V_tip)./lamb(i);
fs = 2 * f_dmax;

[S1,F1,T1,P1] = spectrogram(ephi_i, wlen, hop, nfft, fs); 

F1_doppler = F1 - fs/2; % Doppler frequency [Hz]


S1 = S1./max(S1(:));
S1_norm_db = 20*log(abs(fftshift(S1',2)));
figure(i + 5);

imagesc(F1_doppler*10^(-3), T1, S1_norm_db);

xlabel('Doppler frequency [kHz]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Time[s]', 'FontSize', 12, 'FontWeight', 'bold');
title(['Micro Doppler Specrtum [Db and normalized] N_{fft} = ', num2str(nfft), ' overlap = ', num2str(hop), 'fs = ', num2str(fs)], 'FontSize', 12, 'FontWeight', 'bold');




end
