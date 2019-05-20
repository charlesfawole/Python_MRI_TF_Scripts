function [] = Val_TF_plot_new_304()

% mT_tip = [7.26	6.47	3.39	5.03	3.39	1.49	5.37	4.85	2.80];
% mT_ring = [8.60	7.82	4.06	5.68	3.82	1.68	5.98	5.41	2.95];

mT_tip = [8.33	7.11	3.67	5.48	3.69	1.82	6.31	5.67	3.22];
mT_ring = [9.57	8.72	4.30	6.19	4.12	1.92	6.68	6.26	3.28];

% function rlt = Val_TF(p_coil,ratio, p_gel, p_place, p_TF, p_lead, dy_lead, grid_TF,polate_mode)

%% change these value
p_coil = 'II';
ratio = [1,1];
p_gel = '0';
dx_lead = 0;
dy_lead = 0;
dz_lead = 0;
TFscale_mode = 'mm';


%% temperature tip
cT_tip = Val_TF_new_304(p_coil, ratio, p_gel, 1, 'T', 'tip', dx_lead,dy_lead,dz_lead, 0.2, 'extrap');
%[~,idx] = max(mT_tip);
%cT_tip = cT_tip./cT_tip(idx)*mT_tip(idx);
cT_tip = TFscale(cT_tip, mT_tip, TFscale_mode); % 'relative'
mmax = max([cT_tip mT_tip]);
TF_plot(cT_tip./mmax, mT_tip./mmax,'tip','temperature');
% TF_plot(cT_tip, mT_tip,'tip','temperature');
disp('304 temperature tip')
disp(mT_tip);
disp(cT_tip);


%% temperature ring
cT_ring = Val_TF_new_304(p_coil, ratio, p_gel, 1, 'T', 'ring', dx_lead,dy_lead,dz_lead, 0.2, 'extrap');
% [~,idx] = max(mT_ring);
% cT_ring = cT_ring./cT_ring(idx)*mT_ring(idx);
cT_ring = TFscale(cT_ring, mT_ring, TFscale_mode);
mmax = max([cT_ring mT_ring]);
TF_plot(cT_ring./mmax, mT_ring./mmax,'ring','temperature');
% TF_plot(cT_ring, mT_ring,'ring','temperature');
disp('304 temperature ring')
disp(mT_ring);
disp(cT_ring);

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
title(['Calculated and measured ', t1,' ', t2, ' results (normalized)']);
legend({'Calculated','Measured'})
xlabel('Trajectory Index')
figure
plot(meas,cal,'o'); 
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
xlabel('Measured value')
ylabel('Calculated value')

function cN = TFscale(cN,mN,mode)
switch mode
    case 'abs'
        cN = sum(cN.*mN)/sum(cN.^2)*cN;
    case 'rel'
        cN = sum(cN./mN)/sum((cN./mN).^2)*cN;
    case 'max'
        [~,idx] = max(mN);
        cN = cN./cN(idx)*mN(idx);
    case 'mm' % minimize maximum error
        cN_1 = min(cN./mN);
        cN_2 = max(cN./mN);
        cN = 2/(cN_1+cN_2)*cN;
    otherwise
        [~,idx] = max(mN);
        cN = cN./cN(idx)*mN(idx);
end
