function [] = RegulateTF(name_lead,p_TF,p_lead)

%% change these value
p_coil = 'II';
ratio = [1,1];
p_gel = '0';
p_place = 1;
dy_lead = 0;
polate_mode = 'extrap';
TFscale_mode = 'abs';

%% normalization
dT_simu = 7.97*2.93717/2; 
%Thi is because we use 2W/kg to calculate the dT, while the E field is
%calculated with quadrature excitation with 2.93717 average SAR.
E_rms_simu = 159.7098/sqrt(2);

switch name_lead
    case '303'
        switch p_TF
            case 'T'
                switch p_lead
                    case 'tip'
                        ref_scale = 12.8;
                    case 'ring'
                        ref_scale = 12.59;
                end
            case 'V'
                switch p_lead
                    case 'tip'
                        ref_scale = 127.25;
                    case 'ring'
                        ref_scale = 128.359;
                end
        end
    case '304'
        switch p_TF
            case 'T'
                switch p_lead
                    case 'tip'
                        ref_scale = 13.8;
                    case 'ring'
                        ref_scale = 12.61;
                end
            case 'V'
                switch p_lead
                    case 'tip'
                        ref_scale = 127.25;
                    case 'ring'
                        ref_scale = 128.359;
                end
        end
end
switch p_TF
    case 'T'
        disp('Make sure when study Thermal transfer function, rod temperature should be used as reference.');
        E_scale = sqrt(ref_scale/dT_simu);
    case 'V'
        disp('Make sure when study Thermal transfer function, rod temperature should be used as reference.');
        E_scale = ref_scale/E_rms_simu;
end
    
%%
switch name_lead
    case '303'
        mT_tip = [5.7, 4, 2, 1.83, 1.97, 1.5, 5.34, 3.06, 1.6];
        mT_ring = [5.54, 2.77, 1.6, 1.74, 1.66, 1.6, 5.78, 3.42, 1.42];
        
        mV_tip = [3.55, 3.10, 2.9, 5.33, 3.28, 1.63, 6.29, 4.66, 3.41]./9.39;
        mV_ring = [3.07, 2.77, 1.78, 7.16, 6.44, 4.80, 6.89, 6.23, 4.94]./9.39;
    case '304'
        mT_tip = [9.13, 6.9, 5.63, 7.93, 6.95, 4.46, 14.95, 13.31, 7];
        mT_ring = [7.95, 6.95, 4.16, 5.19, 4.81, 3.2, 12.95, 9.7, 5.59];
        
        mV_tip = [5.98, 3.25, 1.86, 5.79, 4.99, 3.26, 7.23, 6.97, 4.91]./9.39;
        mV_ring = [9.38, 8.65, 7, 6.81, 6.08, 4.61, 9.62, 8.54, 5.29]./9.39;
    otherwise
        disp('The lead is not supported.')
        return;
end

[TF_0,Etan_0,cN] = Val_TF(p_coil, ratio, p_gel, p_place, p_TF, p_lead, dy_lead, 0.2, polate_mode,name_lead);
switch p_TF
    case 'T'
        switch p_lead
            case 'tip'
                mN = mT_tip;
            case 'ring'
                mN = mT_ring;
        end
        p_TF_1 = 'RF heating'
    case 'V'
        switch p_lead
            case 'tip'
                mN = mV_tip;
            case 'ring'
                mN = mV_ring;
        end
        p_TF_1 = 'RF voltage'
end

[r,cN] = TFscale(cN, mN, TFscale_mode);


switch p_TF
    case 'V'
        TF_1 = TF_0.y*r/E_scale;
        disp(abs(TF_1.'));
        disp(angle(TF_1.')*180/pi)
        %disp((abs(Etan_0.y.*E_scale).')/sqrt(2));
        disp(abs(sum(TF_1.*Etan_0.y.*E_scale)));
        %disp(cN);
    case 'T'
        TF_1 = TF_0.y*sqrt(r)/E_scale;
        disp('whole body SAR')
        disp(E_scale^2*2.93717);
        disp('Transfer function')
        disp(abs(TF_1.'));
        disp(angle(TF_1.')*180/pi);
        disp('RMS value of Etan')
        disp((abs(Etan_0.y.*E_scale).')/sqrt(2));
        disp('Magnitude and phase of Electric Field')
        disp((abs(Etan_0.y.*E_scale).'));
        disp((angle(Etan_0.y.*E_scale).')*180/pi);
        %disp(E_scale);
        disp('Calculated temperature rise');
        disp((abs(sum(TF_1.*Etan_0.y.*E_scale)))^2);
        disp('Validate temperature rise');
        disp(cN);
end

figure
plot(TF_0.x*100,abs(TF_1),'b');
grid on
title(['Magnitude of ', p_TF_1, ' transfer function for ', p_lead]);
xlabel('Distance to the tip [cm]');
ylabel('Magnitude of Transfer function');

figure
plot(TF_0.x*100,angle(TF_1)*180/pi,'b');
grid on
title(['Phase of ', p_TF_1, ' transfer function for ',p_lead]);
xlabel('Distance to the tip [cm]');
ylabel('Phase [Degree]');
% legend(p_lead);




function [r,cN] = TFscale(cN,mN,mode)
switch mode
    case 'abs'
        r = sum(cN.*mN)/sum(cN.^2);
    case 'rel'
        r = sum(cN./mN)/sum((cN./mN).^2);
    case 'max'
        [~,idx] = max(mN);
        r = mN(idx)/cN(idx);
    otherwise
        [~,idx] = max(mN);
        cN = cN./cN(idx)*mN(idx);
end
cN = r*cN;


function [TF_0,Etan_0, rlt] = Val_TF(p_coil,ratio, p_gel, p_place, p_TF, p_lead, dy_lead, grid_TF,polate_mode,name_lead)
switch p_coil
    case 'II'
        switch p_gel
            case '+1'
                fid = '3T-II-empty-quad_+10';
            case '0'
                fid = '3T-II-empty-quad';
            case '-1'
                fid = '3T-II-empty-quad_-10';
            otherwise
                disp('The optional parameters are: ''+1'', ''0'' and ''-1''.');
                return;
        end
    case 'III'
        switch p_gel
            case '+1'
                fid = '3T-III-empty-quad_+10';
            case '0'
                fid = '3T-III-empty-quad';
            case '-1'
                fid = '3T-III-empty-quad_-10';
            otherwise
                disp('The optional parameters are: ''+1'', ''0'' and ''-1''.');
                return;
        end
    otherwise
        disp('The optional parameters are: ''II'' and ''III'', please check it.');
        return;
end

Ef = load_ef(fid,ratio);
path = gn_path(p_gel,p_place,p_TF,dy_lead,grid_TF);
grid = search_grid(Ef.coor,path);
Etan = gn_Etan(Ef,grid);


TF = gn_TF(p_lead,p_TF,polate_mode,name_lead);
TF.y =interp1(TF.x,TF.y,(0.5:1:TF.x(end)));
TF.x =(0.5:1:TF.x(end));
TF.x = TF.x/100;

rlt = zeros(1,9);
for ii = 1:1:9
    switch ii
        case 1
            temp = Etan.S1;
        case 2
            temp = Etan.S2;
        case 3
            temp = Etan.S3;
        case 4
            temp = Etan.L1;
        case 5
            temp = Etan.L2;
        case 6
            temp = Etan.L3;
        case 7
            temp = Etan.U1;
        case 8
            temp = Etan.U2;
        case 9
            temp = Etan.U3;
            
    end

    if temp.x(1)>TF.x(1)
        TF.x(1) = [];
        TF.y(1) = [];
    end
    if temp.x(end)<TF.x(end)
        TF.x(end) = [];
        TF.y(end) = [];
    end
    temp.y = interp1(temp.x,temp.y,TF.x);
    temp.x = TF.x;
    if ii ==1
        TF_0 = TF;
        Etan_0 = temp;
    end
    switch p_TF
        case 'V'
            rlt(ii) = abs((temp.y*TF.y.'));
        case 'T'
            rlt(ii) = (abs((temp.y*TF.y.')))^2;
    end
end






%% load the Electric fields
function Ef = load_ef(fid,ratio)
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

Ef.x = ratio(1)*e_x_cos-ratio(2)*1j*e_x_sin;
Ef.y = ratio(1)*e_y_cos-ratio(2)*1j*e_y_sin;
Ef.z = ratio(1)*e_z_cos-ratio(2)*1j*e_z_sin;
Ef.coor.x = x_axis;
Ef.coor.y = y_axis;
Ef.coor.z = z_axis;


%% generate pathways
function path = gn_path(p_gel,p_place,p_TF,dy_lead,grid_TF)
%
if grid_TF<0.2 
    grid_TF = 0.2;
    disp('The grid size of Transfer function is too fine, 0.2 cm grid size is chosen instead.');
end

if grid_TF>1
    grid_TF = 1.0;
    disp('The grid size of Transfer function is too coarse, 1 cm grid size is chosen instead.');
end

%% straight trajectories
path.S1.z = (0:grid_TF:43);
if path.S1.z(end) < 43
    path.S1.z = [path.S1.z 43];
end
path.S1.z = path.S1.z-30+4;
path.S1.x = -21.*ones(size(path.S1.z))+2;
path.S1.y = zeros(size(path.S1.z));
% path.S1.len = sqrt((path.S1.x(2:end)-path.S1.x(1:end-1)).^2+...
%     (path.S1.y(2:end)-path.S1.y(1:end-1)).^2+...
%     (path.S1.z(2:end)-path.S1.z(1:end-1)).^2);

path.S2.x = path.S1.x+4;
path.S2.y = path.S1.y;
path.S2.z = path.S1.z;

path.S3.x = path.S1.x+8;
path.S3.y = path.S1.y;
path.S3.z = path.S1.z;

%% L-shape trajectories
% z direction part
path.L1.z = (0:grid_TF:27);
if path.L1.z(end) < 27
    path.L1.z = [path.L1.z 27];
end
path.L1.z = path.L1.z-30+4;
path.L1.x = -21.*ones(size(path.L1.z))+2;

% x direction part
path.L1.x = [path.L1.x -19+(grid_TF:grid_TF:16)];
if path.L1.x(end) < -19+16
    path.L1.x = [path.L1.x -19+16];
end
path.L1.y = zeros(size(path.L1.x));
path.L1.z = [path.L1.z path.L1.z(end).*ones(1,length(path.L1.x)-length(path.L1.z))];

path.L2.x = path.L1.x+4;
path.L2.y = path.L1.y;
path.L2.z = path.L1.z;

path.L3.x = path.L1.x+8;
path.L3.y = path.L1.y;
path.L3.z = path.L1.z;

%% U-shape trajectories
% z direction part
path.U1.z = (0:grid_TF:17.5);
if path.U1.z(end) < 17.5
    path.U1.z = [path.U1.z 17.5];
end
path.U1.z = path.U1.z-30+4;
path.U1.x = -21.*ones(size(path.U1.z))+2;
% x direction part
path.U1.x = [path.U1.x path.U1.x(end)+(grid_TF:grid_TF:8)];
if path.U1.x(end) < path.U1.x(1)+8
    path.U1.x = [path.U1.x path.U1.x(1)+8];
end
path.U1.z = [path.U1.z path.U1.z(end).*ones(1,length(path.U1.x)-length(path.U1.z))];
% -z direction
path.U1.z = [path.U1.z path.U1.z(end)-(grid_TF:grid_TF:17.5)];
if path.U1.z(end) < path.U1.z(1)
    path.U1.z = [path.U1.z path.U1.z(1)];
end
path.U1.x = [path.U1.x path.U1.x(end).*ones(1,length(path.U1.z)-length(path.U1.x))];
path.U1.y = zeros(size(path.U1.x));
% 
path.U2.x = path.U1.x+4;
path.U2.y = path.U1.y;
path.U2.z = path.U1.z;

path.U3.x = path.U1.x+8;
path.U3.y = path.U1.y;
path.U3.z = path.U1.z;

switch p_place
    case 1
        switch p_TF
            case 'V'
                path.S1.z = -1.*path.S1.z;
                path.S2.z = -1.*path.S2.z;
                path.S3.z = -1.*path.S3.z;
                path.L1.z = -1.*path.L1.z;
                path.L2.z = -1.*path.L2.z;
                path.L3.z = -1.*path.L3.z;
                path.U1.z = -1.*path.U1.z;
                path.U2.z = -1.*path.U2.z;
                path.U3.z = -1.*path.U3.z;
                disp('Trajectories for voltage TF are generated.')
            case 'T'
                path.S1.x = -1.*path.S1.x;
                path.S2.x = -1.*path.S2.x;
                path.S3.x = -1.*path.S3.x;
                path.L1.x = -1.*path.L1.x;
                path.L2.x = -1.*path.L2.x;
                path.L3.x = -1.*path.L3.x;
                path.U1.x = -1.*path.U1.x;
                path.U2.x = -1.*path.U2.x;
                path.U3.x = -1.*path.U3.x;
                disp('Trajectories for temperature TF are generated.')
            otherwise
                disp('TF is not available.');
                path = [];
                return
        end
    case 2
        disp('We will not do this now.')
        return;

    otherwise
        disp('Measurment setting is not available.');
        path = [];
        return;
end

switch p_gel
    case '+1';
        d_y = 1;
    case '-1'
        d_y = -1;
    case '0'
        d_y = 0;
    otherwise
        disp('Gel is not available.')
        path = [];
        return
end

path.S1.y = path.S1.y+d_y+dy_lead;
path.S2.y = path.S2.y+d_y+dy_lead;
path.S3.y = path.S3.y+d_y+dy_lead;
path.L1.y = path.L1.y+d_y+dy_lead;
path.L2.y = path.L2.y+d_y+dy_lead;
path.L3.y = path.L3.y+d_y+dy_lead;
path.U1.y = path.U1.y+d_y+dy_lead;
path.U2.y = path.U2.y+d_y+dy_lead;
path.U3.y = path.U3.y+d_y+dy_lead;

%% change the unit
path.S1.x = path.S1.x/100; path.S1.y = path.S1.y/100; path.S1.z = path.S1.z/100;
path.S2.x = path.S2.x/100; path.S2.y = path.S2.y/100; path.S2.z = path.S2.z/100;
path.S3.x = path.S3.x/100; path.S3.y = path.S3.y/100; path.S3.z = path.S3.z/100;

path.L1.x = path.L1.x/100; path.L1.y = path.L1.y/100; path.L1.z = path.L1.z/100;
path.L2.x = path.L2.x/100; path.L2.y = path.L2.y/100; path.L2.z = path.L2.z/100;
path.L3.x = path.L3.x/100; path.L3.y = path.L3.y/100; path.L3.z = path.L3.z/100;

path.U1.x = path.U1.x/100; path.U1.y = path.U1.y/100; path.U1.z = path.U1.z/100;
path.U2.x = path.U2.x/100; path.U2.y = path.U2.y/100; path.U2.z = path.U2.z/100;
path.U3.x = path.U3.x/100; path.U3.y = path.U3.y/100; path.U3.z = path.U3.z/100;


%% mapping the trajectories into the grid
function grid = search_grid(coor,path)
for ii=1:1:9
    switch ii
        case 1
            traj = path.S1;
        case 2
            traj = path.S2;
        case 3
            traj = path.S3;
        case 4
            traj = path.L1;
        case 5
            traj = path.L2;
        case 6
            traj = path.L3;
        case 7
            traj = path.U1;
        case 8
            traj = path.U2;
        case 9
            traj = path.U3;
    end
    
    idx_grid.x = zeros(size(traj.x));
    idx_grid.y = zeros(size(traj.y));
    idx_grid.z = zeros(size(traj.z));
    for jj = 1:1:length(traj.x)
        idx_grid.x(jj) = find(coor.x>traj.x(jj),1)-1;
        idx_grid.y(jj) = find(coor.y>traj.y(jj),1)-1;
        idx_grid.z(jj) = find(coor.z>traj.z(jj),1)-1;
        if (idx_grid.x(jj)<1) ||(idx_grid.y(jj)<1)||(idx_grid.z(jj)<1)
            disp('Error, the line is out of the gel.')
        end
    end

    kk = 2;
    while kk <= length(idx_grid.x)
        if (idx_grid.x(kk)==idx_grid.x(kk-1))&&(idx_grid.y(kk)==idx_grid.y(kk-1))&&(idx_grid.z(kk)==idx_grid.z(kk-1))
            idx_grid.x(kk) = [];
            idx_grid.y(kk) = [];
            idx_grid.z(kk) = [];
            kk = kk-1;
        end
        kk = kk+1;
    end
    
    switch ii
        case 1
            grid.S1 = idx_grid;
        case 2
            grid.S2 = idx_grid;
        case 3
            grid.S3 = idx_grid;
        case 4
            grid.L1 = idx_grid;
        case 5
            grid.L2 = idx_grid;
        case 6
            grid.L3 = idx_grid;
        case 7
            grid.U1 = idx_grid;
        case 8
            grid.U2 = idx_grid;
        case 9
            grid.U3 = idx_grid;
    end
    
end

%% calculate the Etan
function Etan = gn_Etan(Ef,grid)
for ii = 1:1:9
    switch ii
        case 1
            traj = grid.S1;
        case 2
            traj = grid.S2;
        case 3
            traj = grid.S3;
        case 4
            traj = grid.L1;
        case 5
            traj = grid.L2;
        case 6
            traj = grid.L3;
        case 7
            traj = grid.U1;
        case 8
            traj = grid.U2;
        case 9
            traj = grid.U3;
    end
    x_temp = (Ef.coor.x(traj.x(2:end))-Ef.coor.x(traj.x(1:end-1)));
    y_temp = (Ef.coor.y(traj.y(2:end))-Ef.coor.y(traj.y(1:end-1)));
    z_temp = (Ef.coor.z(traj.z(2:end))-Ef.coor.z(traj.z(1:end-1)));

    x_temp = x_temp.';
    y_temp = y_temp.';
    z_temp = z_temp.';
    
    temp.y =zeros(1,length(x_temp));
    for mm = 1:1:length(traj.x)-1
        dE = [(Ef.x(traj.x(mm+1),traj.y(mm+1),traj.z(mm+1))+Ef.x(traj.x(mm),traj.y(mm),traj.z(mm)))/2
            (Ef.y(traj.x(mm+1),traj.y(mm+1),traj.z(mm+1))+Ef.y(traj.x(mm),traj.y(mm),traj.z(mm)))/2
            (Ef.z(traj.x(mm+1),traj.y(mm+1),traj.z(mm+1))+Ef.z(traj.x(mm),traj.y(mm),traj.z(mm)))/2];
        k = [x_temp(mm),y_temp(mm),z_temp(mm)]./norm([x_temp(mm),y_temp(mm),z_temp(mm)]);
        temp.y(mm) = k*dE;
    end
    
    temp.len = sqrt(x_temp.^2+y_temp.^2+z_temp.^2);
    temp.y = -1.*temp.y(end:-1:1);
    temp.len = temp.len(end:-1:1);
    temp.x = zeros(size(temp.len));
    for ll = 1:1:length(temp.len)
        temp.x(ll) = sum(temp.len(1:1:ll))-temp.len(ll)/2;
    end
    
    switch ii
        case 1
            Etan.S1 = temp;
        case 2
            Etan.S2 = temp;
        case 3
            Etan.S3 = temp;
        case 4
            Etan.L1 = temp;
        case 5
            Etan.L2 = temp;
        case 6
            Etan.L3 = temp;
        case 7
            Etan.U1 = temp;
        case 8
            Etan.U2 = temp;
        case 9
            Etan.U3 = temp;
    end
    
end


%% Generate the Transfer function
function TF = gn_TF(p_lead,p_TF,polate_mode,name_lead)
switch name_lead
    case '303'
        tip2can_T.mag = [0.003458706	0.003459965	0.003223906	0.00299043	0.002762996	0.002551659	0.002391445	0.002289874	0.002265505	0.002269294	0.002352774	0.002448921	0.002572346	0.002721334	0.002841541	0.002913031	0.002954078	0.002986931	0.002980736	0.002903541	0.0028028	0.002665147	0.002477125	0.002261039	0.002000623	0.001712037	0.001396124	0.001068585	0.000804446	0.000679058	0.000710687	0.000900895	0.001224829	0.001470485	0.001699247	0.001922875	0.002128727	0.002267556	0.002435788	0.002512673	0.00256869	0.002564544	0.002535896];
        tip2can_T.phase = [-160.4837494	-160.3305817	-167.3965149	-175.5026398	175.6331024	165.034729	153.3736877	141.896286	130.2715454	117.6055756	106.2945023	96.89645386	86.74334717	78.22200012	72.2793045	67.38743591	62.67749786	58.39102173	54.39930725	50.3474617	47.4747467	44.7316246	41.80550003	38.90817261	35.43113709	30.72391319	23.7647438	11.8422184	-7.882658005	-34.92086411	-62.94311905	-86.63362885	-101.5245895	-108.6764069	-113.1907425	-116.4441147	-118.4347458	-119.4732208	-117.7775421	-115.5310974	-117.9836884	-118.720314	-120.2564774];
        tip2can_T.x = (0:1:length(tip2can_T.mag)-1);
        
        tip2can_V.mag = [0.000465157	0.000464447	0.00046325	0.000459826	0.00045878	0.000453861	0.000446615	0.0004452	0.00044663	0.000447813	0.000451163	0.000453922	0.000454248	0.000452418	0.000444347	0.000436367	0.00042428	0.000408601	0.00038978	0.000371052	0.000345603	0.000320107	0.00029374	0.000267299	0.000245702	0.000230089	0.000217081	0.00021041	0.000209679	0.000215551	0.000231515	0.000249048	0.000265025	0.000281394	0.000295794	0.000308474	0.000314294	0.000316219	0.00031152	0.000301437	0.000284112];
        tip2can_V.phase = [-67.24126434	-66.75868988	-71.31643677	-76.84049225	-83.6499176	-91.1295929	-99.30986023	-108.1896286	-116.1866455	-124.2097092	-131.5852966	-138.0687866	-144.0347748	-149.9198608	-156.0485077	-162.5553436	-168.3851776	-174.5632629	179.4112091	173.2563934	166.8134308	160.0652618	151.8841858	143.365799	133.4320679	122.2576828	111.0633011	99.1481781	85.60111237	71.58562469	58.96666336	48.66986084	40.0159111	33.11748886	26.85930634	21.43569946	16.71414566	13.0012598	9.256832123	5.586853981	1.789436102];
        tip2can_V.mag = tip2can_V.mag(end:-1:1);
        tip2can_V.phase = tip2can_V.phase(end:-1:1);
        tip2can_V.x = (0:1:length(tip2can_V.mag)-1);
        
        ring2can_T.mag = [0.003672028	0.003659072	0.003527627	0.003257754	0.003027298	0.002782114	0.002574409	0.002406021	0.002312976	0.002284955	0.002320461	0.002385258	0.002486171	0.002591722	0.002707839	0.002792636	0.002845355	0.002880624	0.002871431	0.002811128	0.002707651	0.002555616	0.002372277	0.002191341	0.001977057	0.001760494	0.001545779	0.001325533	0.001126589	0.000940669	0.000810207	0.000777608	0.000898078	0.00118792	0.001440082	0.001708618	0.001924605	0.002102378	0.002268355	0.002345724	0.002389362	0.002420429	0.002424517];
        ring2can_T.phase = [-156.210434	-156.3119049	-160.8782959	-168.8273468	-177.2883148	173.0617065	162.6783752	150.7457886	138.8460236	126.0037842	114.4386902	103.4557037	93.70843506	86.10491943	78.64479828	72.21466064	67.62940979	63.65771103	59.79892349	55.49106598	50.80954742	46.21373367	41.20994568	36.65859985	31.66174126	25.57924652	18.9264431	10.87173271	0.476955056	-14.22394276	-33.96760178	-58.13199997	-83.03560638	-98.03648376	-105.7125778	-111.3560486	-114.8820038	-116.7413177	-118.610611	-120.2361374	-121.7271957	-121.7203217	-120.1791611];
        ring2can_T.x = (0:1:length(ring2can_T.mag)-1);
        
        ring2can_V.mag = [0.000606507	0.000587991	0.000638404	0.000686992	0.000727443	0.000767249	0.000781078	0.000790445	0.000782897	0.000763844	0.000742111	0.000703665	0.000674719	0.00063614	0.000606303	0.000586756	0.00057213	0.000570125	0.000571781	0.000590113	0.000608522	0.000640253	0.000674252	0.000708612	0.000752531	0.00079117	0.000838788	0.000891326	0.000938443	0.000988258	0.00101581	0.001060448	0.001077756	0.00110093	0.001110651	0.001113568];
        ring2can_V.phase = [87.73549652	86.8480072	90.57945251	94.00151062	98.16779327	101.8869781	106.1279144	110.3356247	115.0402527	120.37043	125.8569336	133.2083893	140.2946625	149.271759	159.3757935	168.6584015	178.5112915	-171.6044006	-160.8491974	-150.8719177	-141.3014221	-132.506485	-124.1316681	-115.855545	-108.6164551	-101.2935791	-93.9914856	-86.08394623	-78.8370285	-71.52612305	-64.90921783	-57.59177399	-51.63868713	-45.666996	-38.53374863	-33.46471787];
        ring2can_V.x = (0:1:length(ring2can_V.mag)-1);
    case '304'
        tip2can_T.mag = [403	212	163	156	157	159	159	159	150	154	154	149	147	139	134	127	119	110	100	92	83	73	65	58	51	46	44	49	55	61	71	62	88	96	108	118	128	132	140];
        tip2can_T.phase = [-83	-106	-121	-134	-148	-152	-159	-164	-165	-170	-173	-176	-177	178	175	172	168	166	161	157	153	146	139	128	116	98	81	61	53	42	32	27	23	20	16	14	14	12	13];
        tip2can_T.x = (0:1:length(tip2can_T.mag)-1);
        
        tip2can_V.mag = [0.000552828	0.000551404	0.000543298	0.000527555	0.00050612	0.000482502	0.000457961	0.000434798	0.000418237	0.000395703	0.000373396	0.000351544	0.000328998	0.000308184	0.000290193	0.000274164	0.000260664	0.000247827	0.000238538	0.000234239	0.000234792	0.000238216	0.000248434	0.000262575	0.000270915	0.000284973	0.000297269	0.000308446	0.000319944	0.000331089	0.000340338	0.000348167	0.000352792	0.000358116	0.000361472	0.000358356	0.000348	0.000331361	0.000301929	0.000263696];
        tip2can_V.phase = [-61.08677292	-61.17495346	-64.7173996	-69.00744629	-73.46051788	-78.08712006	-82.79831696	-87.22380066	-91.87045288	-96.43453217	-101.0367966	-106.3781815	-111.5713348	-117.7623291	-124.1322021	-130.0285187	-136.8661041	-145.1569214	-154.5147247	-164.0963745	-173.9319458	177.7183075	169.160141	164.1420898	157.6598969	151.5677948	146.7418213	142.8113098	139.2412415	136.348999	132.6543579	129.8197632	126.9333572	124.9659729	123.5341339	121.7962799	118.5855713	115.1712418	110.5344009	106.8135147];
        tip2can_V.x = (0:1:length(tip2can_V.mag)-1);
        tip2can_V.mag = tip2can_V.mag(end:-1:1);
        tip2can_V.phase = tip2can_V.phase(end:-1:1);
        
        ring2can_T.mag = [128	96	91	89	87	88	88	87	87	86	83	82	77	76	71	65	60	54	52	44	42	34	27	25	22	25	23	27	34	38	43	49	54	59	64	69	76	77];
        ring2can_T.phase = [-112	-122	-144	-148	-150	-153	-155	-160	-164	-166	-169	-175	-176	-178	-179	177	171	168	164	160	150	140	138	119	99	85	74	51	48	39	33	29	25	22	25	18	21	18];
        ring2can_T.x = (0:1:length(ring2can_T.mag)-1);
        
        ring2can_V.mag = [1.04	1.18	1.28	1.35	1.42	1.45	1.48	1.48	1.5	1.52	1.52	1.53	1.53	1.52	1.51	1.49	1.49	1.46	1.42	1.37	1.33	1.28	1.25	1.18	1.18	1.15	1.1	1.06	1.03	1.01	1.01	1.03	1.05	1.12	1.16	1.2	1.3	1.4	1.51	1.59];
        ring2can_V.phase = [-10	-5	-2	1	3	4	6	7	8	10	10	12	13	14	16	18	20	22	25	28	31	36	38	42	45	47	53	60	69	79	88	98	104	115	120	126	132	139	147	149];
        ring2can_V.x = (0:1:length(ring2can_V.mag)-1);
end

switch p_lead
    case 'ring'
        switch p_TF
            case 'T'
                TF = ring2can_T;
            case 'V'
                TF = ring2can_V;
            otherwise
                TF = [];
                disp('The optional parameter is ''T'' and ''V''.');
                return
        end
        
    case 'tip'
        switch p_TF
            case 'T'
                TF = tip2can_T;
            case 'V'
                TF = tip2can_V;
            otherwise
                TF = [];
                disp('The optional parameter is ''T'' and ''V''.');
                return
        end
    otherwise
        TF = [];
        disp('The optional parameter is ''ring'' and ''tip''.');
        return
end

switch polate_mode
    case 'interp'
        label.x = TF.x(end);
        label.y1 = TF.mag(end);
        label.y2 = TF.phase(end);
        TF.x = [TF.x 43];
        TF.mag = [TF.mag 0]; TF.phase = [TF.phase TF.phase(end)];
        TF.mag = interp1(TF.x,TF.mag,(0:1:43));
        TF.phase = interp1(TF.x,TF.phase,(0:1:43));
        TF.x = (0:1:43);

    case 'extrap'
        label.x = TF.x(end);
        label.y1 = TF.mag(end);
        label.y2 = TF.phase(end);
        Test.mag = TF.mag;
        Test.phase=TF.phase;
        switch p_lead
            case 'tip'
                TF.mag = interp1(TF.x,TF.mag,(0:1:43),'linear','extrap');
                TF.phase = interp1(TF.x,TF.phase,(0:1:43),'linear','extrap');
                TF.re = interp1(TF.x,real(Test.mag.*exp(1j*Test.phase./180*pi)),(0:1:43),'linear','extrap');
                TF.im = interp1(TF.x,imag(Test.mag.*exp(1j*Test.phase./180*pi)),(0:1:43),'linear','extrap');
                TF.x = (0:1:43);
            case 'ring'
                TF.mag = interp1(TF.x,TF.mag,(0:1:41),'linear','extrap');
                TF.phase = interp1(TF.x,TF.phase,(0:1:41),'linear','extrap');
                TF.re = interp1(TF.x,real(Test.mag.*exp(1j*Test.phase./180*pi)),(0:1:41),'linear','extrap');
                TF.im = interp1(TF.x,imag(Test.mag.*exp(1j*Test.phase./180*pi)),(0:1:41),'linear','extrap');
                TF.x = (0:1:41);
        end

%%        
        TF.mag = abs(TF.re+1j*TF.im);
        TF.phase = angle(TF.re+1j*TF.im)./pi*180;
    case 'no'
        figure
        subplot(2,1,1)
        plot(ring2can_V.x,ring2can_V.mag)
        subplot(2,1,2)
        plot(ring2can_V.x,ring2can_V.phase)
        disp('The interpolation/exterpolation is not performed, raw TF is used.');
    otherwise 
        disp('The parameter is not available.');
        return;
end

TF.y = TF.mag.*exp(1j*TF.phase*pi/180);





