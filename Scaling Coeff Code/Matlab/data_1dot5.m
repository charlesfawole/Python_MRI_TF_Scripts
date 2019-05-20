clear
close all
rms=0.2179;

%% 106 304 T-C
meas_304_tip_V = [3.94 3.26 2.44 2.38 2.18 1.82 1.56 1.04 0.78 2.2]/2*rms;
calc_304_tip_V = [0.4293    0.3605    0.2775    0.2584    0.2156    0.1651    0.1456    0.1009    0.0667 0.18];
meas_304_tip_V=meas_304_tip_V/max(meas_304_tip_V)/1.1;
calc_304_tip_V=calc_304_tip_V/max(calc_304_tip_V);

%% 106 304 R-C
meas_304_ring_V = [3.52 3.00 2.56 2.40 1.74 1.54 1.36  1.04 0.68 2.3]/2*rms;
calc_304_ring_V = [ 0.3836    0.3214    0.2472    0.2362    0.1966    0.1504    0.1193    0.0840    0.0588 0.16];
meas_304_ring_V=meas_304_ring_V/max(meas_304_ring_V)/1.2;
calc_304_ring_V=calc_304_ring_V/max(calc_304_ring_V);

%% 106 303 T-C
meas_303_tip_V = [4.20 3.62 2.88 2.40 2.02 1.92 1.76 1.34 1.02]/2*rms;
calc_303_tip_V = [0.4577    0.3686    0.2663    0.2880    0.2303    0.1653    0.2045    0.1362    0.0770];
meas_303_tip_V = meas_303_tip_V/max(meas_303_tip_V)/1.1;
calc_303_tip_V = calc_303_tip_V/max(calc_303_tip_V);
 
%% 106 303 R-C
meas_303_ring_V = [5.12 4.10 3.08 3.64 2.68 1.88 2.04 1.24 0.98]/2*rms;
calc_303_ring_V = [0.5579    0.4488    0.3241    0.3604    0.2876    0.2063    0.2550    0.1708    0.0981];
meas_303_ring_V = meas_303_ring_V/max(meas_303_ring_V);
calc_303_ring_V = calc_303_ring_V/max(calc_303_ring_V);

meas = [meas_304_tip_V meas_304_ring_V meas_303_tip_V meas_303_ring_V];
calc = [calc_304_tip_V calc_304_ring_V calc_303_tip_V calc_303_ring_V];

corr = sum(meas.*calc)/sqrt(sum(meas.*meas))/sqrt(sum(calc.*calc));
disp('The correlation number of voltage is: ')
disp(corr)




%% 106 304 heating Tip
meas_304_tip_T = ([13.9 9.9 6.0 9.6 6.6 4.0 3.1 2.7 1.4 5.1]);
calc_304_tip_T =[13.9000    9.4604    5.4570   12.0734    8.2251    4.7744    3.9512    2.8066  1.6100 4.01];
meas_304_tip_T=meas_304_tip_T/max(meas_304_tip_T);
calc_304_tip_T=calc_304_tip_T/max(calc_304_tip_T);

%% 106 304 heating Ring
meas_304_ring_T = ([14.3 10.1 5.1 11.8 7.0 3.4 4.8 3.0 1.7 5.4]);
calc_304_ring_T=[14.3000    9.7806    5.6546   11.0949    7.5772    4.3984    2.8600    2.0784    1.2327 4.29];
meas_304_ring_T=meas_304_ring_T/max(meas_304_ring_T)/1.2;
calc_304_ring_T=calc_304_ring_T/max(calc_304_ring_T);

%% 106 303 heating Tip
meas_303_tip_T = ([9.8 6.3 3.2 6.6 3.5 2.2 5.7 4.3 1.9]);
calc_303_tip_T=[9.8000    6.2541    3.2282    7.8621    4.8568    2.5007    9.2327    4.9900    1.7925];
meas_303_tip_T=meas_303_tip_T/max(meas_303_tip_T);
calc_303_tip_T=calc_303_tip_T/max(calc_303_tip_T);

%% 106 303 heating Ring
meas_303_ring_T = ([9.8 7.5 2.7 5.9 4.1 2.5 6.1 4.0 1.8]);
calc_303_ring_T=[9.8000    6.6192    3.6805    3.9420    2.4415    1.2730    9.7520    5.4664    2.1176];
meas_303_ring_T=meas_303_ring_T/max(meas_303_ring_T)/1.1;
calc_303_ring_T=calc_303_ring_T/max(calc_303_ring_T);

meas = [meas_304_tip_T meas_304_ring_T meas_303_tip_T meas_303_ring_T];
calc = [calc_304_tip_T calc_304_ring_T calc_303_tip_T calc_303_ring_T];

corr = sum(meas.*calc)/sqrt(sum(meas.*meas))/sqrt(sum(calc.*calc));
disp('The correlation number of temperature is: ')
disp(corr)
