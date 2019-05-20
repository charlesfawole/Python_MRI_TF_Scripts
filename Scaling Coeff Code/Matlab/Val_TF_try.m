function rlt = Val_TF_try(p_coil,ratio, p_gel, p_place, p_TF, p_lead, dy_lead, grid_TF,polate_mode)
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
% There are two kinds of possible coordinations
% 1. 
%      Left <--> X-; Right <--> X+;
%      Forward <--> Z-; Backward <--> Z+;
%      Bottom <--> Y-; Surface <--> Y+
% 2.   
%      Left <--> X+; Right <--> X-;
%      Forward <--> Z+; Backward <--> Z-;
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


%% heating routines (IPG end position)
% 1. Stright
%    a. 2 cm to right wall, 4 cm to the Forward wall
%    b. 6 cm to right wall, 4 cm to the Forward wall
%    c. 10 cm to right wall, 4 cm to the Forward wall
% 2. L shaple (27 cm backwards along z, then turn left)
%    a. 2 cm to right wall, 4 cm to the Forward wall
%    b. 6 cm to right wall, 4 cm to the Forward wall
%    c. 10 cm to right wall, 4 cm to the Forward wall
% 3. U shaple (17.5 cm backwards along z, then 8 cm turn right, then 17.5 cm forwards)
%    a. 2 cm to right wall, 4 cm to the Forward wall
%    b. 6 cm to right wall, 4 cm to the Forward wall
%    c. 10 cm to right wall, 4 cm to the Forward wall

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


TF = gn_TF(p_lead,p_TF,polate_mode);
TF.y =interp1(TF.x,TF.y,(0.5:1:TF.x(end)));
TF.x =(0.5:1:TF.x(end));
TF.x = TF.x/100;
% figure
% subplot(2,1,1)
% plot(TF.x,abs(TF.y));
% subplot(2,1,2)
% plot(TF.x,angle(TF.y))

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
%             figure
%             subplot(2,1,1)
%             plot(temp.x,abs(temp.y));
%             subplot(2,1,2)
%             plot(temp.x,angle(temp.y))
%             temp.x.'
%             temp.y.'
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
    

    
    switch p_TF
        case 'V'
            rlt(ii) = abs((temp.y*TF.y.')/100);
        case 'T'
            rlt(ii) = (abs((temp.y*TF.y.')/100))^2;
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
%         switch p_TF
%             case 'V'
%                 path.S1.x = -1.*path.S1.x;
%                 path.S2.x = -1.*path.S2.x;
%                 path.S3.x = -1.*path.S3.x;
%                 path.L1.x = -1.*path.L1.x;
%                 path.L2.x = -1.*path.L2.x;
%                 path.L3.x = -1.*path.L3.x;
%                 path.U1.x = -1.*path.U1.x;
%                 path.U2.x = -1.*path.U2.x;
%                 path.U3.x = -1.*path.U3.x;
%                 path.S1.z = -1.*path.S1.z;
%                 path.S2.z = -1.*path.S2.z;
%                 path.S3.z = -1.*path.S3.z;
%                 path.L1.z = -1.*path.L1.z;
%                 path.L2.z = -1.*path.L2.z;
%                 path.L3.z = -1.*path.L3.z;
%                 path.U1.z = -1.*path.U1.z;
%                 path.U2.z = -1.*path.U2.z;
%                 path.U3.z = -1.*path.U3.z;
%                 disp('Trajectories for voltage TF are generated.')
%             case 'T'
%                 path.S1.z = -1.*path.S1.z;
%                 path.S2.z = -1.*path.S2.z;
%                 path.S3.z = -1.*path.S3.z;
%                 path.L1.z = -1.*path.L1.z;
%                 path.L2.z = -1.*path.L2.z;
%                 path.L3.z = -1.*path.L3.z;
%                 path.U1.z = -1.*path.U1.z;
%                 path.U2.z = -1.*path.U2.z;
%                 path.U3.z = -1.*path.U3.z;
%                 disp('Trajectories for temperature TF are generated.')
%             otherwise
%                 disp('TF is not available.');
%                 path = [];
%                 return
%         end
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
function TF = gn_TF(p_lead,p_TF,polate_mode)

tip2can_T.mag = [403	212	163	156	157	159	159	159	150	154	154	149	147	139	134	127	119	110	100	92	83	73	65	58	51	46	44	49	55	61	71	62	88	96	108	118	128	132	140 150 160 170 180];
tip2can_T.phase = [-83	-106	-121	-134	-148	-152	-159	-164	-165	-170	-173	-176	-177	178	175	172	168	166	161	157	153	146	139	128	116	98	81	61	53	42	32	27	23	20	16	14	14	12	13 13 13 13 13 ];
tip2can_T.x = (0:1:length(tip2can_T.mag)-1);
% figure
% subplot(2,1,1)
% plot(tip2can_T.x,tip2can_T.mag)
% subplot(2,1,2)
% plot(tip2can_T.x,tip2can_T.phase)

tip2can_V.mag = [1.12	1.25	1.38	1.475	1.515	1.552	1.554	1.57	1.575	1.55	1.52	1.5	1.45	1.42	1.37	1.31	1.26	1.225	1.19	1.18	1.16	1.13	1.12	1.1	1.08	1.06	1.03	1.03	1.03	1.03	1.04	1.08	1.16	1.24	1.3	1.4	1.49	1.57	1.66	1.75];
tip2can_V.phase = [-9	-3	1	3.8	6	8	9	11	13	16	17.6	20	23	26	29	33	38	42	43	46	48	52	54	56	60	66	69	74	78	87	96	107	118	122	129	136	142	146	150	155];
tip2can_V.x = (0:1:length(tip2can_V.mag)-1);
% figure
% subplot(2,1,1)
% plot(tip2can_V.x,tip2can_V.mag)
% subplot(2,1,2)
% plot(tip2can_V.x,tip2can_V.phase)

ring2can_T.mag = [128	96	91	89	87	88	88	87	87	86	83	82	77	76	71	65	60	54	52	44	42	34	27	25	22	25	23	27	34	38	43	49	54	59	64	69	76	77 81 86 91 96 101];
ring2can_T.phase = [-112	-122	-144	-148	-150	-153	-155	-160	-164	-166	-169	-175	-176	-178	-179	177	171	168	164	160	150	140	138	119	99	85	74	51	48	39	33	29	25	22	25	18	21	18 18 18 18 18 18];
ring2can_T.x = (0:1:length(ring2can_T.mag)-1);
% figure
% subplot(2,1,1)
% plot(ring2can_T.x,ring2can_T.mag)
% subplot(2,1,2)
% plot(ring2can_T.x,ring2can_T.phase)


ring2can_V.mag = [1.04	1.18	1.28	1.35	1.42	1.45	1.48	1.48	1.5	1.52	1.52	1.53	1.53	1.52	1.51	1.49	1.49	1.46	1.42	1.37	1.33	1.28	1.25	1.18	1.18	1.15	1.1	1.06	1.03	1.01	1.01	1.03	1.05	1.12	1.16	1.2	1.3	1.4	1.51	1.59];
ring2can_V.phase = [-10	-5	-2	1	3	4	6	7	8	10	10	12	13	14	16	18	20	22	25	28	31	36	38	42	45	47	53	60	69	79	88	98	104	115	120	126	132	139	147	149];
ring2can_V.x = (0:1:length(ring2can_V.mag)-1);
% figure
% subplot(2,1,1)
% plot(ring2can_V.x,ring2can_V.mag)
% subplot(2,1,2)
% plot(ring2can_V.x,ring2can_V.phase)

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
    case 'inter'
        TF.x = [TF.x 43];
        TF.mag = [TF.mag 0]; TF.phase = [TF.phase TF.phase(end)];
        
        TF.mag = interp1(TF.x,TF.mag,(0:1:43));
        TF.phase = interp1(TF.x,TF.phase,(0:1:43));
        TF.x = (0:1:43);
    case 'extra'
        disp('The exterpolation function is not completed, raw TF is used.')
        
    otherwise
        disp('The interpolation/exterpolation is not performed, raw TF is used.')
end

TF.y = TF.mag.*exp(1j*TF.phase*pi/180);


        


