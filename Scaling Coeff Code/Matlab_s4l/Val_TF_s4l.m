function [rlt,TF,path,Etan] = Val_TF_s4l(fid,ratio, p_gel, p_place, p_TF, p_lead, dx_lead,dy_lead,dz_lead,grid_TF,polate_mode)
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


% Ef = load_ef(fid_cos,fid-sin,ratio);
Ef = load_ef_quadrature(fid);
%path = gn_path(p_gel,p_place,p_TF,dy_lead,grid_TF);
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

rlt = zeros(1,10);
for ii = 1:1:10
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
        case 10 
            temp = Etan.Z1;
            
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


%% load the Electric fields (cos+sin)
function Ef = load_ef(fid_cos,fid_sin,ratio)
load(fid_cos);
Xn = length(Axis0);
Yn = length(Axis1);
Zn = length(Axis2);

e_x_cos = reshape(Snapshot0(:,1),[Xn-1,Yn-1,Zn-1]);
e_y_cos = reshape(Snapshot0(:,2),[Xn-1,Yn-1,Zn-1]);
e_z_cos = reshape(Snapshot0(:,3),[Xn-1,Yn-1,Zn-1]);

clear Snapshot0
clear Axis0 Axis1 Axis2

load(fid_sin);
Xn = length(Axis0);
Yn = length(Axis1);
Zn = length(Axis2);

e_x_sin = reshape(Snapshot0(:,1),[Xn-1,Yn-1,Zn-1]);
e_y_sin = reshape(Snapshot0(:,2),[Xn-1,Yn-1,Zn-1]);
e_z_sin = reshape(Snapshot0(:,3),[Xn-1,Yn-1,Zn-1]);

clear Snapshot0

if prod(size(e_x_cos)==size(e_x_sin))==0
    disp('The sizes of the I/Q Efield are different ');
    return;
end

Ef.x = ratio(1)*e_x_cos-ratio(2)*1j*e_x_sin;
Ef.y = ratio(1)*e_y_cos-ratio(2)*1j*e_y_sin;
Ef.z = ratio(1)*e_z_cos-ratio(2)*1j*e_z_sin;
Ef.coor.x = Axis0.';
Ef.coor.y = Axis1.';
Ef.coor.z = Axis2.';

%% load the Electric fields (quadrature)
function Ef = load_ef_quadrature(fid)
load(fid);
Xn = length(Axis0);
Yn = length(Axis1);
Zn = length(Axis2);

Ef.x = reshape(Snapshot0(:,1),[Xn-1,Yn-1,Zn-1]);
Ef.y = reshape(Snapshot0(:,2),[Xn-1,Yn-1,Zn-1]);
Ef.z = reshape(Snapshot0(:,3),[Xn-1,Yn-1,Zn-1]);

clear Snapshot0

Ef.coor.x = Axis0.';
Ef.coor.y = Axis1.';
Ef.coor.z = Axis2.';


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

%% ZigZag trajectories
zigzag = 'E:\human model\TF\pathway\zigzag.txt';
point = load(zigzag); 
path.Z1.x = point(:,1);
path.Z1.y = point(:,2);
path.Z1.z = point(:,3);

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
    case '+1'
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
for ii=1:1:10
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
        case 10
            traj = path.Z1;
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
        case 10 
            grid.Z1 = idx_grid;
    end
end

%% calculate the Etan
function Etan = gn_Etan(Ef,grid)
for ii = 1:1:10
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
        case 10
            traj = grid.Z1;
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
        case 10
            Etan.Z1 = temp;
    end
    
end


%% Generate the Transfer function
function TF = gn_TF(p_lead,p_TF,polate_mode)

% tip2can_T.mag = [0.001890695	0.001887353	0.001882946	0.001779964	0.001670671	0.001577413	0.001478687	0.001406374	0.001364914	0.001354627	0.001362846	0.001399665	0.001451431	0.001510291	0.001573102	0.001624102	0.001654625	0.001663901	0.001650212	0.001619536	0.001564689	0.001487706	0.001388782	0.001275009	0.001142199	0.000995205	0.000842854	0.000694599	0.000558049	0.000452926	0.000424414	0.000483117	0.000602241	0.000743253	0.000886136	0.001021436];
% tip2can_T.phase = [-41.2722168	-42.30050659	-41.32633972	-49.17824554	-57.18886566	-66.1889801	-77.33455658	-87.83403778	-100.0836487	-110.3128052	-122.1612167	-132.6765289	-142.2387085	-150.9059143	-158.0373077	-165.5658569	-169.6173706	-173.9907837	-178.1348267	178.1137085	173.5793152	170.1749115	165.4567108	161.6906128	157.0620728	151.3158264	144.8917084	135.1526031	121.8022995	101.3906937	74.54537964	49.71170044	32.10609055	21.86352539	15.8217268	12.56908226];
% tip2can_T.x = (0:1:length(tip2can_T.mag)-1);
tip2can_T.mag = [0.003458706	0.003459965	0.003223906	0.00299043	0.002762996	0.002551659	0.002391445	0.002289874	0.002265505	0.002269294	0.002352774	0.002448921	0.002572346	0.002721334	0.002841541	0.002913031	0.002954078	0.002986931	0.002980736	0.002903541	0.0028028	0.002665147	0.002477125	0.002261039	0.002000623	0.001712037	0.001396124	0.001068585	0.000804446	0.000679058	0.000710687	0.000900895	0.001224829	0.001470485	0.001699247	0.001922875	0.002128727	0.002267556	0.002435788	0.002512673	0.00256869	0.002564544	0.002535896];
tip2can_T.phase = [-160.4837494	-160.3305817	-167.3965149	-175.5026398	175.6331024	165.034729	153.3736877	141.896286	130.2715454	117.6055756	106.2945023	96.89645386	86.74334717	78.22200012	72.2793045	67.38743591	62.67749786	58.39102173	54.39930725	50.3474617	47.4747467	44.7316246	41.80550003	38.90817261	35.43113709	30.72391319	23.7647438	11.8422184	-7.882658005	-34.92086411	-62.94311905	-86.63362885	-101.5245895	-108.6764069	-113.1907425	-116.4441147	-118.4347458	-119.4732208	-117.7775421	-115.5310974	-117.9836884	-118.720314	-120.2564774];
tip2can_T.x = (0:1:length(tip2can_T.mag)-1);

% figure
% subplot(2,1,1)
% plot(tip2can_T.x,tip2can_T.mag)
% subplot(2,1,2)
% plot(tip2can_T.x,tip2can_T.phase)


%% tip to can, voltage
% tip2can_V.mag = [0.000217005	0.000216715	0.000236769	0.000256947	0.000274566	0.000286395	0.00029473	0.000298292	0.000296832	0.000290802	0.000281283	0.000268684	0.00025587	0.000241316	0.000226852	0.000219748	0.000213571	0.000206672	0.000205808	0.0002075	0.000213618	0.000222637	0.000234872	0.000248292	0.000260187	0.000275552	0.000291447	0.000307581	0.000328073	0.00034685	0.000367947	0.00038423	0.000401265	0.00041063	0.000419665	0.000426204];
% tip2can_V.phase = [86.68711853	86.74425507	90.58755493	93.91079712	97.48439026	100.9406357	104.6647568	108.132431	112.8381882	118.0619812	122.9589005	129.4950409	136.3804169	144.3047028	155.0171509	163.3429413	172.691864	-177.8878784	-167.8338623	-157.7466888	-147.9480743	-138.4804688	-128.8030396	-120.6410828	-112.6121597	-106.009697	-99.08558655	-92.74256897	-84.45544434	-75.92462158	-68.10050964	-60.19145203	-52.98839951	-45.41592789	-37.43982315	-30.32876968];
% tip2can_V.x = (0:1:length(tip2can_V.mag)-1);

% with 20 cm cable
% tip2can_V.mag = [0.000403656	0.000405867	0.000404667	0.000401468	0.000399372	0.000395974	0.000386879	0.000379553	0.000366178	0.000349177	0.000333455	0.000314312	0.000294002	0.000276566	0.000254707	0.000237017	0.000218764	0.000203372	0.000191221	0.000182525	0.000176957	0.000178813	0.000184592	0.000192609	0.000208229	0.000224381	0.00024107	0.000253763	0.000267144	0.000277083	0.000283635	0.00028168	0.000277429	0.000271021	0.000255376	0.000237224	0.000215596	0.000197652];
% tip2can_V.phase = [-90.8713608	-100.7506027	-109.466156	-116.6771698	-123.6935501	-130.8478546	-137.5808411	-144.1516418	-151.1986389	-158.0426636	-165.3956299	-172.4607697	-179.3491821	173.5471497	166.1092987	157.7663727	149.5105591	139.961441	129.4434662	117.6845856	105.9786758	93.44757843	80.07138062	67.99630737	57.68231583	47.21017075	39.75191498	33.02243805	27.18820763	22.24571419	18.475214	14.14797497	11.03676891	7.835325718	4.673069954	1.701560497	-1.893894076	-4.309242725];

% Sept 27, 2016
% tip2can_V.mag = [0.000166749	0.000172381	0.000170335	0.000159291	0.000147495	0.000135802	0.000128987	0.000125622	0.000124658	0.000123225	0.000125695	0.000128226	0.000127195	0.00012968	0.000129354	0.000127603	0.000124425	0.000121262	0.00011643	0.000108124	9.71565E-05	8.69276E-05	7.55213E-05	6.11869E-05	4.68245E-05	3.66286E-05	2.80881E-05	2.41648E-05	2.87042E-05	4.23892E-05	5.79458E-05	7.05831E-05	8.52363E-05	9.39885E-05	0.000103329	0.000107257	0.000110424	0.000111268	0.000108509	0.00010253	9.51198E-05	0.000089035];
% tip2can_V.phase = [-88.58987427	-91.03997803	-96.42258453	-100.2030945	-107.8698807	-118.220726	-127.7546158	-139.1629486	-148.1694641	-156.927124	-163.9960327	-172.2375488	-179.9105377	174.2230835	168.6120605	163.6905975	157.4130554	153.289505	149.4626617	145.158371	139.8104553	136.2846375	131.9474487	125.6993179	117.87146	103.8135986	81.07175446	47.31532288	6.492486	-13.49111938	-21.178442	-28.39694977	-32.99246597	-36.93267059	-38.86129379	-41.29965973	-43.23103714	-44.79416656	-47.50757599	-49.02203369	-49.24206924	-47.65087128];
% tip2can_V.x = (0:1:length(tip2can_V.mag)-1);
% tip2can_V.mag = tip2can_V.mag(end:-1:1);
% tip2can_V.phase = tip2can_V.phase(end:-1:1);
% tip2can_V.mag(end) = [];
% tip2can_V.phase(end) = [];
% tip2can_V.x(end) = [];

%Sept 29, 2016
%combined
% tip2can_V.mag = [0.000424593	0.000425345	0.000430252	0.000447	0.000464447	0.00046325	0.000459826	0.00045878	0.000453861	0.000446615	0.0004452	0.00044663	0.000447813	0.000451163	0.000453922	0.000454248	0.000452418	0.000444347	0.000436367	0.00042428	0.000408601	0.00038978	0.000371052	0.000345603	0.000320107	0.00029374	0.000267299	0.000245702	0.000230089	0.000217081	0.00021041	0.000209679	0.000215551	0.000231515	0.000249048	0.000265025	0.000281394	0.000295794	0.000308474	0.000314294	0.000316219	0.00031152	0.000301437	0.000284112];
% tip2can_V.phase = [-63.97989273	-64.23647308	-64.34415436	-67.24126434	-66.75868988	-71.31643677	-76.84049225	-83.6499176	-91.1295929	-99.30986023	-108.1896286	-116.1866455	-124.2097092	-131.5852966	-138.0687866	-144.0347748	-149.9198608	-156.0485077	-162.5553436	-168.3851776	-174.5632629	179.4112091	173.2563934	166.8134308	160.0652618	151.8841858	143.365799	133.4320679	122.2576828	111.0633011	99.1481781	85.60111237	71.58562469	58.96666336	48.66986084	40.0159111	33.11748886	26.85930634	21.43569946	16.71414566	13.0012598	9.256832123	5.586853981	1.789436102];
% half into header
tip2can_V.mag = [0.000465157	0.000464447	0.00046325	0.000459826	0.00045878	0.000453861	0.000446615	0.0004452	0.00044663	0.000447813	0.000451163	0.000453922	0.000454248	0.000452418	0.000444347	0.000436367	0.00042428	0.000408601	0.00038978	0.000371052	0.000345603	0.000320107	0.00029374	0.000267299	0.000245702	0.000230089	0.000217081	0.00021041	0.000209679	0.000215551	0.000231515	0.000249048	0.000265025	0.000281394	0.000295794	0.000308474	0.000314294	0.000316219	0.00031152	0.000301437	0.000284112];
tip2can_V.phase = [-67.24126434	-66.75868988	-71.31643677	-76.84049225	-83.6499176	-91.1295929	-99.30986023	-108.1896286	-116.1866455	-124.2097092	-131.5852966	-138.0687866	-144.0347748	-149.9198608	-156.0485077	-162.5553436	-168.3851776	-174.5632629	179.4112091	173.2563934	166.8134308	160.0652618	151.8841858	143.365799	133.4320679	122.2576828	111.0633011	99.1481781	85.60111237	71.58562469	58.96666336	48.66986084	40.0159111	33.11748886	26.85930634	21.43569946	16.71414566	13.0012598	9.256832123	5.586853981	1.789436102];
tip2can_V.mag = tip2can_V.mag(end:-1:1);
tip2can_V.phase = tip2can_V.phase(end:-1:1);
tip2can_V.x = (0:1:length(tip2can_V.mag)-1);
% tip2can_V.mag(end) = [];
% tip2can_V.phase(end) = [];
% tip2can_V.x(end) = [];



% figure
% subplot(2,1,1)
% plot(tip2can_V.x,tip2can_V.mag)
% subplot(2,1,2)
% plot(tip2can_V.x,tip2can_V.phase)

% ring2can_T.mag = [0.00382857	0.003831548	0.003841128	0.003712909	0.003518547	0.003254615	0.002972054	0.002677047	0.002375964	0.002123134	0.001966157	0.001916678	0.001975761	0.002099482	0.00227695	0.002493287	0.00269513	0.002828331	0.002930388	0.00300214	0.003039413	0.003042819	0.00297786	0.002855058	0.002669596	0.002417084	0.002135751	0.001797518	0.001441817	0.001079813	0.000762072	0.000615067	0.000755632	0.00104455	0.001277148	0.001494995];
% ring2can_T.phase = [-37.76859665	-37.81377411	-37.60815048	-40.15307236	-45.47552872	-51.53127289	-58.50635529	-66.73126221	-77.46324921	-90.47071838	-106.0938263	-122.5702591	-137.5870361	-150.6467133	-161.3736877	-170.3991852	-177.1580353	178.5149841	174.7002106	171.6533508	168.3126068	164.4228668	160.8714142	157.4510956	154.1174927	150.8457336	146.8639832	141.9673615	135.1316833	124.2319336	104.3244705	68.0168457	29.73788834	10.30328178	1.790606022	-3.294094086];
% ring2can_T.x = (0:1:length(ring2can_T.mag)-1);

ring2can_T.mag = [0.003672028	0.003659072	0.003527627	0.003257754	0.003027298	0.002782114	0.002574409	0.002406021	0.002312976	0.002284955	0.002320461	0.002385258	0.002486171	0.002591722	0.002707839	0.002792636	0.002845355	0.002880624	0.002871431	0.002811128	0.002707651	0.002555616	0.002372277	0.002191341	0.001977057	0.001760494	0.001545779	0.001325533	0.001126589	0.000940669	0.000810207	0.000777608	0.000898078	0.00118792	0.001440082	0.001708618	0.001924605	0.002102378	0.002268355	0.002345724	0.002389362	0.002420429	0.002424517];
ring2can_T.phase = [-156.210434	-156.3119049	-160.8782959	-168.8273468	-177.2883148	173.0617065	162.6783752	150.7457886	138.8460236	126.0037842	114.4386902	103.4557037	93.70843506	86.10491943	78.64479828	72.21466064	67.62940979	63.65771103	59.79892349	55.49106598	50.80954742	46.21373367	41.20994568	36.65859985	31.66174126	25.57924652	18.9264431	10.87173271	0.476955056	-14.22394276	-33.96760178	-58.13199997	-83.03560638	-98.03648376	-105.7125778	-111.3560486	-114.8820038	-116.7413177	-118.610611	-120.2361374	-121.7271957	-121.7203217	-120.1791611];
ring2can_T.x = (0:1:length(ring2can_T.mag)-1);
% figure
% subplot(2,1,1)
% plot(ring2can_T.x,ring2can_T.mag)
% subplot(2,1,2)
% plot(ring2can_T.x,ring2can_T.phase)


ring2can_V.mag = [0.000606507	0.000587991	0.000638404	0.000686992	0.000727443	0.000767249	0.000781078	0.000790445	0.000782897	0.000763844	0.000742111	0.000703665	0.000674719	0.00063614	0.000606303	0.000586756	0.00057213	0.000570125	0.000571781	0.000590113	0.000608522	0.000640253	0.000674252	0.000708612	0.000752531	0.00079117	0.000838788	0.000891326	0.000938443	0.000988258	0.00101581	0.001060448	0.001077756	0.00110093	0.001110651	0.001113568];
ring2can_V.phase = [87.73549652	86.8480072	90.57945251	94.00151062	98.16779327	101.8869781	106.1279144	110.3356247	115.0402527	120.37043	125.8569336	133.2083893	140.2946625	149.271759	159.3757935	168.6584015	178.5112915	-171.6044006	-160.8491974	-150.8719177	-141.3014221	-132.506485	-124.1316681	-115.855545	-108.6164551	-101.2935791	-93.9914856	-86.08394623	-78.8370285	-71.52612305	-64.90921783	-57.59177399	-51.63868713	-45.666996	-38.53374863	-33.46471787];
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
        figure
        subplot(2,1,1)
        plot(TF.x,TF.mag)
        hold on
        scatter(label.x,label.y1);
        subplot(2,1,2)
        plot(TF.x,TF.phase);
        hold on
        scatter(label.x,label.y2);
        disp('The interpolation function is used.');
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
        figure
        subplot(2,1,1)
        plot(TF.x,TF.mag,'b')
        hold on
        plot(TF.x,abs(TF.re+1j*TF.im),'r');
        hold on
        plot(label.x,label.y1,'o');
        subplot(2,1,2)
        plot(TF.x,TF.phase,'b');
        hold on
        plot(TF.x,angle(TF.re+1j*TF.im)./pi*180,'r');
        hold on
        plot(label.x,label.y2,'o');
        disp('The exterpolation function is used.');
%        
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

%% ZigZag trajectories
zigzag = 'E:\human model\TF\pathway\zigzag.txt';
point = load(zigzag); 
path.Z1.x = point(:,1);
path.Z1.y = point(:,2);
path.Z1.z = point(:,3);

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
    case '+1'
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
  


