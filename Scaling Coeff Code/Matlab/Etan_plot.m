function Etan = Etan_plot(p_coil,ratio, p_gel, dx_lead,dy_lead, dz_lead, grid_TF)
%% simulations
% We have 6*2 quadratrue simulations in total, as:
% type II coil
%     case 1: Gel [-210,-50,-300], [210,50,300]
%     case 2(+10): Gel [-210,-40,-300],[210,60,300]
%     case 3(-10): Gel [-210,-60,-200],[210,60,300]
% type III coil (rotate 22.5 degree comparing with type II)
%     case 1: Gel [-210,-50,-300], [210,50,300]
%     case 2(+10): Gel [-210,-40,-300],[210,60,300]
%     case 3(-10): Gel [-210,-60,-200],[210,60,300]


%% cordinations
% There are one kinds of possible coordinations
% 1. 
%      Left <--> X-; Right <--> X+;
%      Forward <--> Z-; Backward <--> Z+;
%      Bottom <--> Y-; Surface <--> Y+



%% voltage routines (IPG end position)
% 1. Stright
%    a. 2 cm to left wall, 4 cm to the Forward wall
%    b. 6 cm to left wall, 4 cm to the Forward wall
%    c. 10 cm to left wall, 4 cm to the Forward wall
% 2. L shaple (27 cm backwards along z, then turn right)
%    a. 2 cm to left wall, 4 cm to the Forward wall
%    b. 6 cm to left wall, 4 cm to the Forward wall
%    c. 10 cm to left wall, 4 cm to the Forward wall
% 3. U shaple (17.5 cm backwards along z, then 8 cm turn right, then 17.5 cm forwards)
%    a. 2 cm to left wall, 4 cm to the Forward wall
%    b. 6 cm to left wall, 4 cm to the Forward wall
%    c. 10 cm to left wall, 4 cm to the Forward wall


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
path = gn_path(p_gel,dx_lead,dy_lead,dz_lead,grid_TF);
grid_path = search_grid(Ef.coor,path);
Etan = gn_Etan(Ef,grid_path);

rr = 1/sqrt(2.93717/2);
figure
plot(Etan.S1.x*100,rr*abs(Etan.S1.y),'b');
hold on
plot(Etan.L1.x*100,rr*abs(Etan.L1.y),'r');
plot(Etan.U1.x*100,rr*abs(Etan.U1.y),'g');
grid on
title('Magnitude');
xlabel('Distance to the tip [cm]');
ylabel('Magnitude [V/m]');
legend('Straight path','L shape path','U shape path');

figure
plot(Etan.S1.x*100,angle(Etan.S1.y)*180/pi,'b');
hold on
plot(Etan.L1.x*100,angle(Etan.L1.y)*180/pi,'r');
plot(Etan.U1.x*100,angle(Etan.U1.y)*180/pi,'g');
grid on
title('Phase');
xlabel('Distance to the tip [cm]');
ylabel('Phase [Degree]');
legend('Straight path','L shape path','U shape path');

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
function path = gn_path(p_gel,dx_lead,dy_lead,dz_lead,grid_TF)
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

path.S1.x = path.S1.x+dx_lead;
path.S2.x = path.S2.x+dx_lead;
path.S3.x = path.S3.x+dx_lead;
path.L1.x = path.L1.x+dx_lead;
path.L2.x = path.L2.x+dx_lead;
path.L3.x = path.L3.x+dx_lead;
path.U1.x = path.U1.x+dx_lead;
path.U2.x = path.U2.x+dx_lead;
path.U3.x = path.U3.x+dx_lead;

path.S1.z = path.S1.z+dz_lead;
path.S2.z = path.S2.z+dz_lead;
path.S3.z = path.S3.z+dz_lead;
path.L1.z = path.L1.z+dz_lead;
path.L2.z = path.L2.z+dz_lead;
path.L3.z = path.L3.z+dz_lead;
path.U1.z = path.U1.z+dz_lead;
path.U2.z = path.U2.z+dz_lead;
path.U3.z = path.U3.z+dz_lead;

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



        


