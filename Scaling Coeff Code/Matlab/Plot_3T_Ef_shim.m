function [] = Plot_3T_Ef_shim(fid,y_coor,ratio)
%% The measured results

data = [-19	-25	61.2
-11	-25	68.4
-3	-25	77.8
-19	-15	102.6
-11	-15	73.3
-3	-15	55.9
-19	-5	118
-11	-5	85.6
-3	-5	36
-19	5	114
-11	5	80.1
-3	5	33.94
-19	15	99.5
-11	15	65.4
-3	15	43.7
-19	25	62
-11	25	63.3
-3	25	67.5
19	-25	68.3
11	-25	68.7
3	-25	78.1
19	-15	118.5
11	-15	81.6
3	-15	50.1
19	-5	145.89
11	-5	98.7
3	-5	31
19	5	144
11	5	98.7
3	5	31.1
19	15	122.3
11	15	86.4
3	15	47.6
19	25	69
11	25	64.3
3	25	65.9];

%% load the simulation Electric field and compute the combined Electrical field.
% fid_e = 'Efield_-10mm.mat';
%y_coor = -1;
fid_cos = [fid,'-cos.mat'];
load(fid_cos);
e_x_cos = Real_Part_of_X_Comp_0s + 1j*Imaginary_Part_of_X_Comp_0s;
e_y_cos = Real_Part_of_Y_Comp_0s + 1j*Imaginary_Part_of_Y_Comp_0s;
e_z_cos = Real_Part_of_Z_Comp_0s + 1j*Imaginary_Part_of_Z_Comp_0s;
clear Real_Part_of_X_Comp_0s Imaginary_Part_of_X_Comp_0s
clear Real_Part_of_Y_Comp_0s Imaginary_Part_of_Y_Comp_0s
clear Real_Part_of_Z_Comp_0s Imaginary_Part_of_Z_Comp_0s
clear x_axis y_axis z_axis

fid_sin = [fid,'-sin.mat'];
load(fid_sin);
e_x_sin = Real_Part_of_X_Comp_0s + 1j*Imaginary_Part_of_X_Comp_0s;
e_y_sin = Real_Part_of_Y_Comp_0s + 1j*Imaginary_Part_of_Y_Comp_0s;
e_z_sin = Real_Part_of_Z_Comp_0s + 1j*Imaginary_Part_of_Z_Comp_0s;
clear Real_Part_of_X_Comp_0s Imaginary_Part_of_X_Comp_0s
clear Real_Part_of_Y_Comp_0s Imaginary_Part_of_Y_Comp_0s
clear Real_Part_of_Z_Comp_0s Imaginary_Part_of_Z_Comp_0s

if prod(size(e_x_cos)==size(e_x_sin))==0
    disp('The sizes of the I/Q Efield are different ');
    return;
end

e_x = ratio(1)*e_x_cos-ratio(2)*1j*e_x_sin;
e_y = ratio(1)*e_y_cos-ratio(2)*1j*e_y_sin;
e_z = ratio(1)*e_z_cos-ratio(2)*1j*e_z_sin;

e_rms_all = sqrt(abs(e_x).^2+abs(e_y).^2+abs(e_z).^2)/sqrt(2);


%% load the measured results
x_coor = data(:,1)./100;% cm
z_coor = data(:,2)./100; % cm
y_coor = y_coor./1000; % mm
e_rms = data(:,3);


%% 
e_rms_simu = e_rms.*0;
[~,iy] = min(abs(y_axis-y_coor));
for ii = 1:1:length(e_rms)
    [~,ix] = min(abs(x_axis-x_coor(ii)));
    [~,iz] = min(abs(z_axis-z_coor(ii)));
    e_rms_simu(ii) = e_rms_all(ix,iy,iz);
end

n_mode = 'l';
switch n_mode
    case 'l'
        b1 = e_rms_simu\e_rms;
        e_rms_cal = b1*e_rms_simu;
    case 's'
        e_rms_cal = sqrt(sum(e_rms.^2)/sum(e_rms_simu.^2)).*e_rms_simu;
end

corr = (e_rms_simu.'*e_rms)/sqrt(e_rms_simu.'*e_rms_simu)/sqrt(e_rms.'*e_rms);
err_r = (e_rms_cal-e_rms)./e_rms;
err_abs = (e_rms_cal-e_rms);
disp('correlation is');
disp(corr);
disp('relative error is');
disp(err_r);
disp('absolute error is');
disp(err_abs);

% if strcmp(f_mode,'on')
%     figure
%     scatter(e_rms_cal,e_rms)
%     hold on
%     yl = max(e_rms)-min(e_rms);
%     plot(e_rms_cal,e_rms_cal,'- r')
%     %plot(e_rms_cal,e_rms_cal*1.3,'- m')
%     %plot(e_rms_cal,e_rms_cal*0.7,'- m')
%     xlim([min(e_rms)-0.05*yl,max(e_rms)+0.05*yl]);
%     ylim([min(e_rms)-0.05*yl,max(e_rms)+0.05*yl])
%     
%     xlabel('Normalized Simulation Electric Field Strength (V/m)')
%     ylabel('Measured Electric Field Strength (V/m)')
%     title('Comparison Between Simulations and Experiments of Electric Field Strength')
%     grid on
%     
%     %Rsq1 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2);
%     %Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
%     txt1 = ['correlation coefficient is: ',num2str(corr)];
%     txt2 = ['std of relative error is: ',num2str(std(err_r)*100),'%'];
%     txt3 = ['maximum relative error is: ',num2str(max(abs(err_r))*100),'%'];
%     %txt2 = ['R^2 for slope is: ', num2str(Rsq1)];
%     %txt3 = ['R^2 for slope & Intercept is: ', num2str(Rsq2)];
%     text(min(e_rms),max(e_rms),txt1)
%     text(min(e_rms),max(e_rms)-0.075*yl,txt2)
%     text(min(e_rms),max(e_rms)-0.15*yl,txt3)
% end



e_rms_plane = zeros(length(x_axis),length(z_axis));
e_x_plane = zeros(length(x_axis),length(z_axis));
e_y_plane = zeros(length(x_axis),length(z_axis));
e_z_plane = zeros(length(x_axis),length(z_axis));

for iz = 1:1:length(z_axis)
    e_rms_plane(:,iz) = e_rms_all(:,iy,iz);
    e_x_plane(:,iz) = e_x(:,iy,iz);
    e_y_plane(:,iz) = e_y(:,iy,iz);
    e_z_plane(:,iz) = e_z(:,iy,iz);
end
e_rms_plane = b1*e_rms_plane;
e_x_plane = b1*e_x_plane;
e_y_plane = b1*e_y_plane;
e_z_plane = b1*e_z_plane;


%%
for ix = 1:1:length(x_axis)
    if x_axis(ix)>=-0.210
        x_axis(1:ix) = [];
        e_rms_plane(1:ix,:) = [];
        e_x_plane(1:ix,:) = [];
        e_y_plane(1:ix,:) = [];
        e_z_plane(1:ix,:) = [];
        break;
    end
end
for ix = length(x_axis):-1:1
    if x_axis(ix)<=0.210
        x_axis(ix:end) = [];
        e_rms_plane(ix:end,:) = [];
        e_x_plane(ix:end,:) = [];
        e_y_plane(ix:end,:) = [];
        e_z_plane(ix:end,:) = [];
        break;
    end
end
%%
for iz = 1:1:length(z_axis)
    if z_axis(iz)>=-0.300
        z_axis(1:iz) = [];
        e_rms_plane(:,1:iz) = [];
        e_x_plane(:,1:iz) = [];
        e_y_plane(:,1:iz) = [];
        e_z_plane(:,1:iz) = [];
        break;
    end
end
for iz = length(z_axis):-1:1
    if z_axis(iz)<=0.300
        z_axis(iz:end) = [];
        e_rms_plane(:,iz:end) = [];
        e_x_plane(:,iz:end) = [];
        e_y_plane(:,iz:end) = [];
        e_z_plane(:,iz:end) = [];
        break;
    end
end


figure
surf(z_axis,x_axis,e_rms_plane,'LineStyle','None')
hold on
for ii = 1:1:length(e_rms)
    [~,ix] = min(abs(x_axis-x_coor(ii)));
    [~,iz] = min(abs(z_axis-z_coor(ii)));
    scatter3(z_axis(iz),x_axis(ix),e_rms(ii),'filled')
end
hold off

% check
figure
e_rms_check = sqrt(abs(e_x_plane).^2+abs(e_y_plane).^2+abs(e_z_plane).^2)./sqrt(2);
%e_rms_all = sqrt(abs(e_x).^2+abs(e_y).^2+abs(e_z).^2)/sqrt(2);

disp(max(max(max(abs(e_rms_check-e_rms_plane)))));
surf(z_axis,x_axis,e_rms_check,'LineStyle','None');
hold on
for ii = 1:1:length(e_rms)
    [~,ix] = min(abs(x_axis-x_coor(ii)));
    [~,iz] = min(abs(z_axis-z_coor(ii)));
    scatter3(z_axis(iz),x_axis(ix),e_rms(ii),'filled')
end
hold off
save('Ef_plane.mat','e_x_plane','e_y_plane','e_z_plane','x_axis','z_axis');



