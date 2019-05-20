function [] = Val_TF_plot_304()

mT_tip = [9.13, 6.9, 5.63, 7.93, 6.95, 4.46, 14.95, 13.31, 7];
mT_ring = [7.95, 6.95, 4.16, 5.19, 4.81, 3.2, 12.95, 9.7, 5.59];

mV_tip = [2.39, 0.79, 0.55, 3.17, 2.48, 1.44, 4.04, 3.29, 2.27];
mV_ring = [9.38, 8.65, 7, 6.81, 6.08, 4.61, 9.62, 8.54, 5.29];

% function rlt = Val_TF(p_coil,ratio, p_gel, p_place, p_TF, p_lead, dy_lead, grid_TF,polate_mode)

%% change these value
p_coil = 'II';
ratio = [1,1];
p_gel = '0';
dy_lead = 0;


%% temperature tip
cT_tip = Val_TF(p_coil, ratio, p_gel, 1, 'T', 'tip', dy_lead, 0.2, 'extrap');
[~,idx] = max(mT_tip);
cT_tip = cT_tip./cT_tip(idx)*mT_tip(idx);
TF_plot(cT_tip, mT_tip,'tip','temperature');


%% temperature ring
cT_ring = Val_TF(p_coil, ratio, p_gel, 1, 'T', 'ring', dy_lead, 0.2, 'extrap');
[~,idx] = max(mT_ring);
cT_ring = cT_ring./cT_ring(idx)*mT_ring(idx);
TF_plot(cT_ring, mT_ring,'ring','temperature');

%% voltage tip
cV_tip = Val_TF(p_coil, ratio, p_gel, 1, 'V', 'tip', dy_lead, 0.2, 'extrap');
[~,idx] = max(mV_tip);
cV_tip = cV_tip./cV_tip(idx)*mV_tip(idx);
TF_plot(cV_tip, mV_tip,'tip','voltage');


%% voltage ring
cV_ring = Val_TF(p_coil, ratio, p_gel, 1, 'V', 'ring', dy_lead, 0.2, 'extrap');
[~,idx] = max(mV_ring);
cV_ring = cV_ring./cV_ring(idx)*mV_ring(idx);
TF_plot(cV_ring, mV_ring,'ring','voltage');

function [] = TF_plot(cal,meas,t1,t2)
% figure
% plot(cal./meas);
% grid on;
% names = {'S1'; 'S2'; 'S3'; 'L1'; 'L2'; 'L3'; 'U1'; 'U2'; 'U3'};
% set(gca,'xtick',(1:9),'xticklabel',names);
% title(['Calculated/Measured ',t1,' ',t2,' results results']);


figure
b = bar([cal.', meas.']);
b(1).FaceColor = 'blue';
b(2).FaceColor = 'red';
grid on;
names = {'S1'; 'S2'; 'S3'; 'L1'; 'L2'; 'L3'; 'U1'; 'U2'; 'U3'};
set(gca,'xtick',(1:9),'xticklabel',names);
title(['Calculated and measured ', t1,' ', t2, ' results']);
legend({'Calculated','Measured'})
