function [] = Val_Zig_304_plot()

mT_tip = [9.13, 6.9, 5.63, 7.93, 6.95, 4.46, 14.95, 13.31, 7]; mT_tip_z = 3.5;
mT_ring = [7.95, 6.95, 4.16, 5.19, 4.81, 3.2, 12.95, 9.7, 5.59]; mT_ring_z = 2.9;

% mV_tip = [5.98, 3.25, 0.84, 5.79, 4.99, 3.26, 7.23, 6.97, 4.91];
% mV_tip = [2.39, 0.79, 0.55, 3.17, 2.48, 1.44, 4.04, 3.29, 2.27]; mV_tip_z = 1.06;
mV_tip = [5.98, 3.25, 1.86, 5.79, 4.99, 3.26, 7.23, 6.97, 4.91]; mV_tip_z = 2.85;
%mV_tip = [4.68, 3.59, 1.86, 3.17, 2.48, 1.44, 4.04, 3.29, 2.27]; mV_tip_z = 1.06;
mV_ring = [9.38, 8.65, 7, 6.81, 6.08, 4.61, 9.62, 8.54, 5.29]; mV_ring_z = 5.55;

% function rlt = Val_TF(p_coil,ratio, p_gel, p_place, p_TF, p_lead, dy_lead, grid_TF,polate_mode)

%% change these value
p_coil = 'II';
ratio = [1,1];
p_gel = '0';
dx_lead = 0;
dy_lead = 0;
dz_lead = 0;
TFscale_mode = 'abs';


%% temperature tip
cT_tip = Val_Zig_304(p_coil, ratio, p_gel, 1, 'T', 'tip', dx_lead, dy_lead, dz_lead, 0.2, 'extrap');
%[~,idx] = max(mT_tip);
%cT_tip = cT_tip./cT_tip(idx)*mT_tip(idx);
[rr,~] = TFscale(cT_tip(1:9), mT_tip, TFscale_mode); % 'relative'
cT_tip = cT_tip.*rr;
mT_tip = [mT_tip mT_tip_z];
mmax = max([cT_tip mT_tip]);
TF_plot(cT_tip/mmax, mT_tip/mmax,'tip','temperature');
% TF_plot(cT_tip, mT_tip,'tip','temperature');
disp('304 temperature tip')
disp(mT_tip);
disp(cT_tip);

%% temperature ring
cT_ring = Val_Zig_304(p_coil, ratio, p_gel, 1, 'T', 'ring', dx_lead, dy_lead, dz_lead, 0.2, 'extrap');
% [~,idx] = max(mT_ring);
% cT_ring = cT_ring./cT_ring(idx)*mT_ring(idx);
[rr,~] = TFscale(cT_ring(1:9), mT_ring, TFscale_mode);
cT_ring = cT_ring.*rr;
mT_ring = [mT_ring mT_ring_z];
mmax = max([cT_ring mT_ring]);
TF_plot(cT_ring./mmax, mT_ring./mmax,'ring','temperature');
% TF_plot(cT_ring, mT_ring,'ring','temperature');
disp('304 temperature ring')
disp(mT_ring);
disp(cT_ring);


%% voltage tip
cV_tip = Val_Zig_304(p_coil, ratio, p_gel, 1, 'V', 'tip', dx_lead, dy_lead, dz_lead, 0.2, 'extrap');
% [~,idx] = max(mV_tip);
% cV_tip = cV_tip./cV_tip(idx)*mV_tip(idx);
[rr,~] = TFscale(cV_tip(1:9), mV_tip, TFscale_mode);
cV_tip = cV_tip.*rr;
mV_tip = [mV_tip mV_tip_z];
mmax = max([cV_tip mV_tip]);
TF_plot(cV_tip./mmax, mV_tip./mmax,'tip','voltage');
% TF_plot(cV_tip, mV_tip,'tip','voltage');
disp('304 voltage tip')
disp(mV_tip);
disp(cV_tip);

% voltage ring
cV_ring = Val_Zig_304(p_coil, ratio, p_gel, 1, 'V', 'ring', dx_lead, dy_lead, dz_lead, 0.2, 'extrap');
% [~,idx] = max(mV_ring);
% cV_ring = cV_ring./cV_ring(idx)*mV_ring(idx);
[rr,~] = TFscale(cV_ring(1:9), mV_ring, TFscale_mode);
cV_ring = cV_ring.*rr;
mV_ring = [mV_ring mV_ring_z];
mmax = max([cV_ring mV_ring]);
TF_plot(cV_ring./mmax, mV_ring./mmax,'ring','voltage');
% TF_plot(cV_ring, mV_ring,'ring','voltage');
disp('304 voltage ring')
disp(mV_ring);
disp(cV_ring);


function [] = TF_plot(cal,meas,t1,t2)
% max
% mmax = max([cal meas]);
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
names = {'S1'; 'S2'; 'S3'; 'L1'; 'L2'; 'L3'; 'U1'; 'U2'; 'U3';'Z'};
set(gca,'xtick',(1:10),'xticklabel',names);
title(['Calculated and measured ', t1,' ', t2, ' results (normalized)']);
legend({'Calculated','Measured'})
xlabel('Trajectory Index')
figure
plot(cal,meas,'o'); 
hold on
axis([0,max(max(cal),max(meas)),0, max(max(cal),max(meas))]);
plot([0,max([max(cal),max(meas)])],[0,max([max(cal),max(meas)])])
hold on
switch t2
    case 'temperature'
        plot([0,max([max(cal),max(meas)])],[0,(1-2*0.2242)*max([max(cal),max(meas)])],'g');
        plot([0,max([max(cal),max(meas)])],[0,(1+2*0.2242)*max([max(cal),max(meas)])],'g');
    case 'voltage'
        plot([0,max([max(cal),max(meas)])],[0,(1-2*0.23)*max([max(cal),max(meas)])],'g');
        plot([0,max([max(cal),max(meas)])],[0,(1+2*0.23)*max([max(cal),max(meas)])],'g');
end

title(['Correlation plot of ', t1,' ', t2, ' results (normalized)']);
xlabel('Calculated value')
ylabel('Measured value')

function [rr,cN] = TFscale(cN,mN,mode)
switch mode
    case 'abs'
        rr = sum(cN.*mN)/sum(cN.^2);
    case 'rel'
        rr = sum(cN./mN)/sum((cN./mN).^2);
    case 'max'
        [~,idx] = max(mN);
        rr = 1/cN(idx)*mN(idx);
    otherwise
        [~,idx] = max(mN);
        rr = 1/cN(idx)*mN(idx);
end
cN = rr.*cN;