function rlt = Val_TF_304(p_coil,ratio, p_gel, p_place, p_TF, p_lead, dx_lead, dy_lead, dz_lead, grid_TF,polate_mode)
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
path = gn_path2(p_gel,p_place,p_TF,dx_lead,dy_lead,dz_lead,grid_TF);
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
%     figure
%     subplot(2,1,1)
%     plot(temp.x,abs(temp.y));
%     title(['the E_{tan} of ',num2str(ii)]);
%     subplot(2,1,2)
%     plot(temp.x,angle(temp.y).*180/pi);
    
end


%% Generate the Transfer function
function TF = gn_TF(p_lead,p_TF,polate_mode)

tip2can_T.mag = [403	212	163	156	157	159	159	159	150	154	154	149	147	139	134	127	119	110	100	92	83	73	65	58	51	46	44	49	55	61	71	62	88	96	108	118	128	132	140];
tip2can_T.phase = [-83	-106	-121	-134	-148	-152	-159	-164	-165	-170	-173	-176	-177	178	175	172	168	166	161	157	153	146	139	128	116	98	81	61	53	42	32	27	23	20	16	14	14	12	13];
tip2can_T.x = (0:1:length(tip2can_T.mag)-1);
% figure
% subplot(2,1,1)
% plot(tip2can_T.x,tip2can_T.mag)
% subplot(2,1,2)
% plot(tip2can_T.x,tip2can_T.phase)

% previous 
% tip2can_V.mag = [1.12	1.25	1.38	1.475	1.515	1.552	1.554	1.57	1.575	1.55	1.52	1.5	1.45	1.42	1.37	1.31	1.26	1.225	1.19	1.18	1.16	1.13	1.12	1.1	1.08	1.06	1.03	1.03	1.03	1.03	1.04	1.08	1.16	1.24	1.3	1.4	1.49	1.57	1.66	1.75];
% tip2can_V.phase = [-9	-3	1	3.8	6	8	9	11	13	16	17.6	20	23	26	29	33	38	42	43	46	48	52	54	56	60	66	69	74	78	87	96	107	118	122	129	136	142	146	150	155];
% tip2can_V.x = (0:1:length(tip2can_V.mag)-1);
%% with 20 cm cable
% tip2can_V.mag = [0.000604273	0.000606233	0.000571444	0.000536617	0.000496702	0.000462068	0.000431241	0.000398199	0.000366746	0.000340724	0.000320061	0.00030229	0.000287191	0.00027885	0.000275111	0.000275628	0.000287604	0.000298878	0.000321474	0.000337925	0.000359305	0.000381583	0.000399684	0.00041713	0.00043076	0.000442073	0.000453905	0.000458709	0.000461949	0.000462417	0.000461052	0.000454882	0.000444382	0.000424208	0.000383524	0.000320611	0.000269884	0.000233528];
% tip2can_V.phase = [-71.55413055	-71.65146637	-76.00460052	-81.20170593	-85.92933655	-91.25571442	-96.3677063	-102.3378983	-108.6606216	-115.6484604	-122.8180847	-131.3025055	-140.4827728	-150.5174255	-159.427948	-168.9516754	-178.8888397	172.8226013	165.0442505	158.8935699	152.7797852	148.2072449	144.3424683	140.6292725	137.596817	135.2241058	133.0065613	130.9339905	129.1394348	127.4095001	126.1519547	124.1736298	122.6995316	120.4507141	118.0440445	113.4227905	109.0964432	106.8304749];

%% head removed Sept 27, 2016
% tip2can_V.mag = [0.000259657	0.000259472	0.000249602	0.000241014	0.00023259	0.000227022	0.000217993	0.000208961	0.00020345	0.00019641	0.000185602	0.00017619	0.000169344	0.000160252	0.000152411	0.000139206	0.000131412	0.00012209	0.00011015	9.88835E-05	8.80784E-05	7.98076E-05	6.72316E-05	5.78046E-05	4.69678E-05	3.57092E-05	2.69853E-05	1.88093E-05	1.61143E-05	1.95791E-05	2.30357E-05	3.18312E-05	3.44832E-05	3.96416E-05	4.21442E-05	4.85501E-05	4.93794E-05	4.78465E-05	5.27443E-05	5.42317E-05	0.000053301	4.95173E-05];
% tip2can_V.phase = [-104.1095963	-104.6341858	-108.1249619	-112.1911926	-115.8786011	-120.0298233	-121.731514	-126.0800629	-129.3321533	-132.8068542	-134.5424652	-139.6641998	-142.4179077	-144.7249298	-147.4250946	-149.5997162	-152.1529083	-153.589859	-155.946228	-157.1149597	-161.0305634	-162.2870331	-166.1314697	-169.385788	-175.1361542	179.8621063	168.8313141	145.3509674	110.4082031	80.48234558	65.52493286	59.56293869	47.95347595	44.69169617	37.76108551	35.16250229	34.84330368	31.01413155	34.4210434	32.04130554	27.25432396	27.4231987];
%% Sept 29, 2016
tip2can_V.mag = [0.000552828	0.000551404	0.000543298	0.000527555	0.00050612	0.000482502	0.000457961	0.000434798	0.000418237	0.000395703	0.000373396	0.000351544	0.000328998	0.000308184	0.000290193	0.000274164	0.000260664	0.000247827	0.000238538	0.000234239	0.000234792	0.000238216	0.000248434	0.000262575	0.000270915	0.000284973	0.000297269	0.000308446	0.000319944	0.000331089	0.000340338	0.000348167	0.000352792	0.000358116	0.000361472	0.000358356	0.000348	0.000331361	0.000301929	0.000263696];
tip2can_V.phase = [-61.08677292	-61.17495346	-64.7173996	-69.00744629	-73.46051788	-78.08712006	-82.79831696	-87.22380066	-91.87045288	-96.43453217	-101.0367966	-106.3781815	-111.5713348	-117.7623291	-124.1322021	-130.0285187	-136.8661041	-145.1569214	-154.5147247	-164.0963745	-173.9319458	177.7183075	169.160141	164.1420898	157.6598969	151.5677948	146.7418213	142.8113098	139.2412415	136.348999	132.6543579	129.8197632	126.9333572	124.9659729	123.5341339	121.7962799	118.5855713	115.1712418	110.5344009	106.8135147];

%% Oct 5, 2016
% tip2can_V.mag = [0.000884265	0.000882946	0.000860206	0.000817768	0.00080082	0.00077783	0.000757182	0.000718541	0.000679419	0.000628477	0.000588464	0.000556908	0.000523067	0.000493342	0.000457708	0.00041904	0.000391785	0.000364177	0.000343385	0.000327998	0.000321102	0.000318851	0.000328241	0.00033263	0.000345887	0.000362437	0.000378815	0.00039774	0.000413789	0.000433693	0.000452204	0.000468668	0.00047977	0.00048988	0.000500443	0.000503502	0.000503957	0.000502768	0.000495003	0.000488009	0.00046803	0.00042708	0.000377917];
% tip2can_V.phase = [142.2742157	142.2670593	139.748703	130.0687866	125.9664993	121.9570007	118.4460678	113.8110886	108.9723358	103.6888733	98.48791504	94.37786102	90.17624664	84.78939056	79.68769073	73.51351929	66.64718628	58.41381454	49.22366714	39.4425621	29.63975716	20.69767952	13.12988949	3.772838593	-3.376145601	-9.98637867	-15.55441093	-21.36186218	-26.25545883	-29.965065	-33.65948486	-36.76945877	-39.53726959	-41.8779335	-43.88668442	-46.16165924	-47.62870407	-50.05149078	-51.89680862	-53.25894165	-55.12284851	-58.76196289	-62.03815842];

tip2can_V.x = (0:1:length(tip2can_V.mag)-1);
tip2can_V.mag = tip2can_V.mag(end:-1:1);
tip2can_V.phase = tip2can_V.phase(end:-1:1);

% figure
% subplot(2,1,1)
% plot(tip2can_V.x,tip2can_V.mag)
% subplot(2,1,2)
% plot(tip2can_V.x,tip2can_V.phase)

ring2can_T.mag = [128	96	91	89	87	88	88	87	87	86	83	82	77	76	71	65	60	54	52	44	42	34	27	25	22	25	23	27	34	38	43	49	54	59	64	69	76	77];
ring2can_T.phase = [-112	-122	-144	-148	-150	-153	-155	-160	-164	-166	-169	-175	-176	-178	-179	177	171	168	164	160	150	140	138	119	99	85	74	51	48	39	33	29	25	22	25	18	21	18];
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
    case 'interp'
        label.x = TF.x(end);
        label.y1 = TF.mag(end);
        label.y2 = TF.phase(end);
        TF.x = [TF.x 43];
        TF.mag = [TF.mag 0]; TF.phase = [TF.phase TF.phase(end)];
        TF.mag = interp1(TF.x,TF.mag,(0:1:43));
        TF.phase = interp1(TF.x,TF.phase,(0:1:43));
        TF.x = (0:1:43);
%         figure
%         subplot(2,1,1)
%         plot(TF.x,TF.mag)
%         hold on
%         scatter(label.x,label.y1);
%         subplot(2,1,2)
%         plot(TF.x,TF.phase);
%         hold on
%         scatter(label.x,label.y2);
%         disp('The interpolation function is used.');
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
%         figure
%         subplot(2,1,1)
%         plot(TF.x,TF.mag,'b')
%         hold on
%         plot(TF.x,abs(TF.re+1j*TF.im),'r')
%         scatter(label.x,label.y1);
%         subplot(2,1,2)
%         plot(TF.x,TF.phase,'b');
%         hold on
%         plot(TF.x,angle(TF.re+1j*TF.im)./pi*180,'r')
%         scatter(label.x,label.y2);
%         disp('The exterpolation function is used.');
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

%% generate pathways
function path = gn_path2(p_gel,p_place,p_TF,dx_lead,dy_lead,dz_lead,grid_TF)
%
if grid_TF<0.2 
    grid_TF = 0.2;
    disp('The grid size of Transfer function is too fine, 0.2 cm grid size is chosen instead.');
end

if grid_TF>1
    grid_TF = 1.0;
    disp('The grid size of Transfer function is too coarse, 1 cm grid size is chosen instead.');
end

if abs(dx_lead)>2 ||abs(dy_lead)>2 ||abs(dz_lead)>2
    disp('The number is too large, please check it. The unit is CM')
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
%% dx
path.S1.x = path.S1.x+dx_lead;
path.S2.x = path.S2.x+dx_lead;
path.S3.x = path.S3.x+dx_lead;
path.L1.x = path.L1.x+dx_lead;
path.L2.x = path.L2.x+dx_lead;
path.L3.x = path.L3.x+dx_lead;
path.U1.x = path.U1.x+dx_lead;
path.U2.x = path.U2.x+dx_lead;
path.U3.x = path.U3.x+dx_lead;
%% dy
path.S1.y = path.S1.y+d_y+dy_lead;
path.S2.y = path.S2.y+d_y+dy_lead;
path.S3.y = path.S3.y+d_y+dy_lead;
path.L1.y = path.L1.y+d_y+dy_lead;
path.L2.y = path.L2.y+d_y+dy_lead;
path.L3.y = path.L3.y+d_y+dy_lead;
path.U1.y = path.U1.y+d_y+dy_lead;
path.U2.y = path.U2.y+d_y+dy_lead;
path.U3.y = path.U3.y+d_y+dy_lead;
% dz
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
        


