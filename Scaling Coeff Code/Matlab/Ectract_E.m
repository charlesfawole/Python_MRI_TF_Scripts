function [] = Ectract_E(p_coil, ratio, p_gel, axis,pos)

pos = pos/100;
%%
switch axis
    case 'x'
        if (pos(1)>5||pos(1)<-5)||(pos(2)>30||pos(2)<-30)
            disp('The size of the gel is 42*10*60 (x*y*z cm^3)');
            disp('y/z is in range of [-5,5]/[-30,30] cm');
            return
        end
    case 'y'
        if (pos(1)>21||pos(1)<-21)||(pos(2)>30||pos(2)<-30)
            disp('The size of the gel is 42*10*60 (x*y*z cm^3)');
            disp('x/z is in range of [-21,21]/[-30,30] cm');
            return
        end
    case 'z'
        if (pos(1)>21||pos(1)<-21)||(pos(2)>5||pos(2)<-5)
            disp('The size of the gel is 42*10*60 (x*y*z cm^3)');
            disp('y/z is in range of [-5,5]/[-30,30] cm');
            return
        end
    otherwise
        disp('The optional parameters are ''x'', ''y'', and ''z''');

end


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

switch axis
    case 'x'
        idx_y = find(Ef.coor.y>pos(1),1)-1;
        idx_z = find(Ef.coor.z>pos(2),1)-1;
        Ef_plot.x = Ef.x(:,idx_y,idx_z);
        Ef_plot.y = Ef.y(:,idx_y,idx_z);
        Ef_plot.z = Ef.z(:,idx_y,idx_z);
        ax = Ef.coor.x.';
    case 'y'
        idx_x = find(Ef.coor.x>pos(1),1)-1;
        idx_z = find(Ef.coor.z>pos(2),1)-1;
        Ef_plot.x = Ef.x(idx_x,:,idx_z);
        Ef_plot.y = Ef.y(idx_x,:,idx_z);
        Ef_plot.z = Ef.z(idx_x,:,idx_z);
        ax = Ef.coor.y.';
    case 'z'
        idx_x = find(Ef.coor.x>pos(1),1)-1;
        idx_y = find(Ef.coor.y>pos(2),1)-1;
        aa = size(Ef.x);
        Ef_plot.x = zeros(aa(3),1);
        Ef_plot.y = zeros(aa(3),1);
        Ef_plot.z = zeros(aa(3),1);
        for ii = 1:1:aa(3)
            Ef_plot.x(ii) = Ef.x(idx_x,idx_y,ii);
            Ef_plot.y(ii) = Ef.y(idx_x,idx_y,ii);
            Ef_plot.z(ii) = Ef.z(idx_x,idx_y,ii);
        end
        ax = Ef.coor.z.';
        
end

figure
title('x component of Efield along the line');
subplot(2,1,1)
plot(ax,abs(Ef_plot.x))
grid on
hold on
switch axis
    case 'x'
        scatter([-0.21,0.21] ,[0, 0]);
    case 'y'
        scatter([-0.05,0.05] ,[0, 0]);
    case 'z'
        scatter([-0.30,0.30] ,[0, 0]);
        scatter([-0.225,0.225] ,[0, 0]);
end
subplot(2,1,2)

plot(ax,unwrap(angle(Ef_plot.x)))
grid on
hold on
switch axis
    case 'x'
        scatter([-0.21,0.21] ,[0, 0]);
    case 'y'
        scatter([-0.05,0.05] ,[0, 0]);
    case 'z'
        scatter([-0.30,0.30] ,[0, 0]);
        scatter([-0.225,0.225] ,[0, 0]);
end


figure
title('y component of Efield along the line');
subplot(2,1,1)
plot(ax,abs(Ef_plot.y))
grid on
hold on
switch axis
    case 'x'
        scatter([-0.21,0.21] ,[0, 0]);
    case 'y'
        scatter([-0.05,0.05] ,[0, 0]);
    case 'z'
        scatter([-0.30,0.30] ,[0, 0]);
        scatter([-0.225,0.225] ,[0, 0]);
end
subplot(2,1,2)
plot(ax,unwrap(angle(Ef_plot.y)))
grid on
hold on
switch axis
    case 'x'
        scatter([-0.21,0.21] ,[0, 0]);
    case 'y'
        scatter([-0.05,0.05] ,[0, 0]);
    case 'z'
        scatter([-0.30,0.30] ,[0, 0]);
        scatter([-0.225,0.225] ,[0, 0],'r');
end

figure
title('z component of Efield along the line');
subplot(2,1,1)
plot(ax,abs(Ef_plot.z))
grid on
hold on
switch axis
    case 'x'
        scatter([-0.21,0.21] ,[0, 0]);
    case 'y'
        scatter([-0.05,0.05] ,[0, 0]);
    case 'z'
        scatter([-0.30,0.30] ,[0, 0]);
        scatter([-0.225,0.225] ,[0, 0],'r');
end
subplot(2,1,2)
plot(ax,unwrap(angle(Ef_plot.z)))
hold on
switch axis
    case 'x'
        scatter([-0.21,0.21] ,[0, 0]);
    case 'y'
        scatter([-0.05,0.05] ,[0, 0]);
    case 'z'
        scatter([-0.225,0.225] ,[0, 0],'r');
end

figure
title('Total Efield along the line');
plot(ax,sqrt(abs(Ef_plot.x).^2+abs(Ef_plot.y).^2+abs(Ef_plot.z).^2));
hold on
switch axis
    case 'x'
        scatter([-0.21,0.21] ,[0, 0]);
    case 'y'
        scatter([-0.05,0.05] ,[0, 0]);
    case 'z'
        scatter([-0.225,0.225] ,[0, 0],'r');
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



