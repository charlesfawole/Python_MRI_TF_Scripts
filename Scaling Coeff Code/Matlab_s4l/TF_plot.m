function [] = TF_plot(cal,meas,t1,t2)
% figure
% plot(cal./meas);
% grid on;
% names = {'S1'; 'S2'; 'S3'; 'L1'; 'L2'; 'L3'; 'U1'; 'U2'; 'U3'};
% set(gca,'xtick',(1:9),'xticklabel',names);
% title(['Calculated/Measured ',t1,' ',t2,' results results']);
figure
b = bar([cal.', meas.']);
%b(1).FaceColor = 'blue';
%b(2).FaceColor = 'red';
grid on;
names = {'S1'; 'S2'; 'S3'; 'L1'; 'L2'; 'L3'; 'U1'; 'U2'; 'U3';'Z1'};
set(gca,'xtick',(1:10),'xticklabel',names);
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
        plot([0,max([max(cal),max(meas)])],[0,(1-2*0.2242)*max([max(cal),max(meas)])],'g'); %% one sigma is 22.42% for 3T temperature & 23.81% for 1.5T temperature
        plot([0,max([max(cal),max(meas)])],[0,(1+2*0.2242)*max([max(cal),max(meas)])],'g');
    case 'voltage'
        plot([0,max([max(cal),max(meas)])],[0,(1-2*0.23)*max([max(cal),max(meas)])],'g');% one simga is 23.01% for 3T voltage and 14.66% for 1.5T voltage
        plot([0,max([max(cal),max(meas)])],[0,(1+2*0.23)*max([max(cal),max(meas)])],'g');
end
title(['Correlation plot of ', t1,' ', t2, ' results (normalized)']);
xlabel('Measured value')
ylabel('Calculated value')
end