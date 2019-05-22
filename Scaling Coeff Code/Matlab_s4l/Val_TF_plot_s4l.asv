%function [] = Val_TF_plot_303_s4l()
clear
mT_tip = [5.7, 4, 2, 1.83, 1.97, 1.5, 5.34, 3.06, 1.6, 4]; %% the last number for zigzag test only
mT_ring = [5.54, 2.77, 1.6, 1.74, 1.66, 1.6, 5.78, 3.42, 1.42];

%mV_tip = [3.55, 3.10, 2.9, 5.33, 3.28, 1.63, 6.29, 4.66, 3.41];
mV_tip = [3.55, 3.10, 2.9, 5.33, 3.28, 4.31, 6.29, 4.66, 3.41];
mV_ring = [3.07, 2.77, 1.78, 7.16, 6.44, 4.80, 6.89, 6.23, 4.94];

% function rlt = Val_TF(p_coil,ratio, p_gel, p_place, p_TF, p_lead, dy_lead, grid_TF,polate_mode)

%% change these value
fid = '3T-coil-II-ZigZag_quadrature_SEMX14.mat';
% fid_cos = '3T-coil-II-ZigZag_cos_SEMX14';
% fid_sin = '3T-coil-II-ZigZag_sin_SEMX14';
ratio = [1,1];
p_gel = '0';
dx_lead = 0;
dy_lead = 0;
dz_lead = 0;
TFscale_mode = 'max';
p_place = 1;
p_TF = 'T';
p_lead = 'tip';
grid_TF = 0.2;
polate_mode = 'extrap';

%% temperature tip
[cT_tip,TF,path,Etan] = Val_TF_s4l(fid, ratio, p_gel, p_place, p_TF, p_lead, dx_lead,dy_lead,dz_lead,grid_TF,polate_mode);
%[~,idx] = max(mT_tip);
%cT_tip = cT_tip./cT_tip(idx)*mT_tip(idx);
[cT_tip,TF,TF_result] = TFscale(p_TF,cT_tip, mT_tip, TFscale_mode,TF); % 'relative'
mmax = max([cT_tip mT_tip]);
TF_plot(cT_tip./mmax, mT_tip./mmax,'tip','temperature');
% TF_plot(cT_tip, mT_tip,'tip','temperature');
disp('temperature tip')
disp(mT_tip);
disp(cT_tip);


% %% temperature ring
% cT_ring = Val_TF_303(p_coil, ratio, p_gel, 1, 'T', 'ring', dx_lead,dy_lead,dz_lead, 0.2, 'extrap');
% % [~,idx] = max(mT_ring);
% % cT_ring = cT_ring./cT_ring(idx)*mT_ring(idx);
% cT_ring = TFscale(cT_ring, mT_ring, TFscale_mode);
% mmax = max([cT_ring mT_ring]);
% TF_plot(cT_ring./mmax, mT_ring./mmax,'ring','temperature');
% % TF_plot(cT_ring, mT_ring,'ring','temperature');
% disp('303 temperature ring')
% disp(mT_ring);
% disp(cT_ring);
% 
% %% voltage tip
% cV_tip = Val_TF_303(p_coil, ratio, p_gel, 1, 'V', 'tip', dx_lead,dy_lead,dz_lead, 0.2, 'extrap');
% % [~,idx] = max(mV_tip);
% % cV_tip = cV_tip./cV_tip(idx)*mV_tip(idx);
% cV_tip = TFscale(cV_tip, mV_tip, TFscale_mode);
% mmax = max([cV_tip mV_tip]);
% TF_plot(cV_tip./mmax, mV_tip./mmax,'tip','voltage');
% % TF_plot(cV_tip, mV_tip,'tip','voltage');
% disp('303 voltage tip')
% disp(mV_tip);
% disp(cV_tip);
% 
% 
% %% voltage ring
% cV_ring = Val_TF_303(p_coil, ratio, p_gel, 1, 'V', 'ring', dx_lead,dy_lead,dz_lead, 0.2, 'extrap');
% % [~,idx] = max(mV_ring);
% % cV_ring = cV_ring./cV_ring(idx)*mV_ring(idx);
% cV_ring = TFscale(cV_ring, mV_ring, TFscale_mode);
% mmax = max([cV_ring mV_ring]);
% TF_plot(cV_ring./mmax, mV_ring./mmax,'ring','voltage');
% % TF_plot(cV_ring, mV_ring,'ring','voltage');
% disp('303 voltage ring')
% disp(mV_ring);
% disp(cV_ring);



