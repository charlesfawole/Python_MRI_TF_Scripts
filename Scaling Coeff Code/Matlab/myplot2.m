function [] = myplot2()
%% voltage

mode = 'abs';
rms = 0.2179;
% % 304
% mV_tip = [3.94 3.26 2.44 2.38 2.18 1.82 1.56 1.04 0.78]/2*rms;
% cV_tip = [0.4293 0.3605 0.2775 0.2584 0.2156 0.1651 0.1456 0.1009 0.0667];
% cV_tip = TFscale(cV_tip,mV_tip,mode)
% 
% mV_ring = [3.52 3.00 2.56 2.40 1.74 1.54 1.36 1.04 0.68]/2*rms;
% cV_ring = [0.3836 0.3214 0.2472 0.2362 0.1966 0.1504 0.1193 0.0840 0.0588];
% cV_ring = TFscale(cV_ring,mV_ring,mode)

% % 303
% 
% mV_tip = [4.20 3.62 2.88 2.40 2.02 1.92 1.76 1.34 1.02]/2*rms;
% cV_tip = [0.4577    0.3686    0.2663    0.2880    0.2303    0.1653    0.2045    0.1362    0.0770];
% % cV_tip = TFscale(cV_tip,mV_tip,mode)

% mV_ring = [5.12 4.10 3.08 3.64 2.68 1.88 2.04 1.24 0.98]/2*rms;
% cV_ring = [0.5579    0.4488    0.3241    0.3604    0.2876    0.2063    0.2550    0.1708    0.0981];
% cV_ring = TFscale(cV_ring,mV_ring,mode)

% heating
% 303
% mT_tip = [9.8 6.3 3.2 6.6 3.5 2.2 5.7 4.3 1.9];
% cT_tip = [9.8000    6.2541    3.2282    7.8621    4.8568    2.5007    9.2327    4.9900    1.7925];
% cT_tip = TFscale(cT_tip,mT_tip,mode)
% 
mT_ring = [9.8 7.5 2.7 5.9 4.1 2.5 6.1 4.0 1.8];
cT_ring = [9.8000    6.6192    3.6805    3.9420    2.4415    1.2730    9.7520    5.4664    2.1176];
cT_ring = TFscale(cT_ring,mT_ring,mode);

%304
% mT_tip = [13.9 9.9 6.0 9.6 6.6 4.0 3.1 2.7 1.4];
% cT_tip = [13.9000    9.4604    5.4570   12.0734    8.2251    4.7744    3.9512    2.8066  1.6100];
% cT_tip = TFscale(cT_tip,mT_tip,mode);
% 
% mT_ring = [14.3 10.1 5.1 11.8 7.0 3.4 4.8 3.0 1.7];
% cT_ring = [14.3000    9.7806    5.6546   11.0949    7.5772    4.3984    2.8600    2.0784    1.2327];
% cT_ring = TFscale(cT_ring,mT_ring,mode);

%
% mmax = max([cV_tip mV_tip]);
% TF_plot(cV_tip./mmax, mV_tip./mmax,'tip','voltage');
% 
% mmax = max([cV_ring mV_ring]);
% TF_plot(cV_ring./mmax, mV_ring./mmax,'ring','voltage');

% mmax = max([cT_tip mT_tip]);
% TF_plot(cT_tip./mmax, mT_tip./mmax,'tip','temperature');

mmax = max([cT_ring mT_ring]);
TF_plot(cT_ring./mmax, mT_ring./mmax,'ring','temperature');

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
names = {'S1'; 'S2'; 'S3'; 'L1'; 'L2'; 'L3'; 'U1'; 'U2'; 'U3'};
set(gca,'xtick',(1:9),'xticklabel',names);
title(['Calculated and measured ', t1,' ', t2, ' results (normalized)']);
legend({'Calculated','Measured'})
figure
plot(cal,meas,'o'); 
hold on
axis([0,max(max(cal),max(meas)),0, max(max(cal),max(meas))]);
plot((0:max([max(cal),max(meas)])),(0:max([max(cal),max(meas)])))
hold on
plot([0,max([max(cal),max(meas)])],[0,0.5*max([max(cal),max(meas)])],'g');
plot([0,max([max(cal),max(meas)])],[0,1.5*max([max(cal),max(meas)])],'g');
title(['Correlation plot of ', t1,' ', t2, ' results']);

function cN = TFscale(cN,mN,mode)
switch mode
    case 'abs'
        cN = sum(cN.*mN)/sum(cN.^2)*cN;
    case 'rel'
        cN = sum(cN./mN)/sum((cN./mN).^2)*cN;
    case 'max'
        [~,idx] = max(mN);
        cN = cN./cN(idx)*mN(idx);
    otherwise
        [~,idx] = max(mN);
        cN = cN./cN(idx)*mN(idx);
end