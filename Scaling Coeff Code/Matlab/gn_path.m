%% generate pathways
function path = gn_path(p_gel,p_place,p_TF,grid_TF)
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

path.S2.x =path.S1.x+4;
path.S2.y = path.S1.y;
path.S2.z = path.S1.z;

path.S3.x =path.S1.x+8;
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

path.L2.x =path.L1.x+4;
path.L2.y = path.L1.y;
path.L2.z = path.L1.z;

path.L3.x =path.L1.x+8;
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
path.U2.x =path.U1.x+4;
path.U2.y = path.U1.y;
path.U2.z = path.U1.z;

path.U3.x =path.U1.x+8;
path.U3.y = path.U1.y;
path.U3.z = path.U1.z;

switch p_place
    case 1
        switch p_TF
            case 'V'
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
                disp('Trajectories for voltage TF are generated.')
            otherwise
                disp('TF is not available.');
                path = [];
                return
        end
    case 2
        switch p_TF
            case 'V'
                path.S1.x = -1.*path.S1.x;
                path.S2.x = -1.*path.S2.x;
                path.S3.x = -1.*path.S3.x;
                path.L1.x = -1.*path.L1.x;
                path.L2.x = -1.*path.L2.x;
                path.L3.x = -1.*path.L3.x;
                path.U1.x = -1.*path.U1.x;
                path.U2.x = -1.*path.U2.x;
                path.U3.x = -1.*path.U3.x;
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
            otherwise
                disp('TF is not available.');
                path = [];
                return
        end
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

path.S1.y = path.S1.y+d_y;
path.S2.y = path.S2.y+d_y;
path.S3.y = path.S3.y+d_y;
path.L1.y = path.L1.y+d_y;
path.L2.y = path.L2.y+d_y;
path.L3.y = path.L3.y+d_y;
path.U1.y = path.U1.y+d_y;
path.U2.y = path.U2.y+d_y;
path.U3.y = path.U3.y+d_y;

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