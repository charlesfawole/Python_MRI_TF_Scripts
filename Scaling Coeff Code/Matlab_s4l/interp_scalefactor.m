TF.position = flipud(xlsread('104 to 106.xlsx','3T','A2:A37'));
TF.mag0 = xlsread('104 to 106.xlsx','3T','E2:E37');
TF.phase0 = xlsread('104 to 106.xlsx','3T','I2:I37');
TF.complex = TF.mag0.*exp(j*((TF.phase0/180)*pi));

TF.y = interp1(TF.position,TF.complex,(1:1:43),'linear','extrap');  %TF.y is complex
TF.x = 0.01:0.01:0.43;
TF.mag = abs(TF.y);
TF.phase = angle(TF.y);
TF.phase = TF.phase./pi.*180;
result_file = 'interpTF104_3T.csv';
TF_result = [TF.mag.',TF.phase.'];
csvwrite(result_file,TF_result);
%xlswrite(result_file,TF_result,strcat('A1',':','B',num2str(length(TF.phase))));

figure
subplot(2,1,1)
plot(TF.x*100,TF.mag,'b')
hold on
plot(TF.position,TF.mag0,'r')
hold off
title('extrapolated transfer function magnitude')
subplot(2,1,2)
plot(TF.x*100,TF.phase,'b');
hold on
plot(TF.position,TF.phase0,'r')
hold off
title('extrapolated transfer function phase')


p_TF = 'T';
mea = []; % measured temperature or voltage
Edata = 'quadrature_SEMX14_S2.txt';
Efield = load(Edata);
Etan.y = Efield(:,4)+1j*Efield(:,5);
Etan.y = -1*Etan.y(end:-1:1);
Etan.d = [0;(sqrt(sum(diff(Efield(:,1:3)).^2,2)))];
Etan.x = [];
Etan.x(1) = 0;
for i = 1:1:length(Etan.d)-1
    Etan.x(i+1) = Etan.x(i)+Etan.d(i+1);
end
leadlength = length(TF.y)/100;
TF.x = 0.01:0.01:leadlength;
Etan.y=interp1(Etan.x,Etan.y,TF.x,'linear','extrap');

switch p_TF
    case 'V'
        scaling_factor = abs((Etan.y*TF.y.')/100)
    case 'T'
        scaling_factor = (abs((Etan.y*TF.y.')/100))^2
end

deltaT_304_M106 = 6;


scaledTF = TF.y*(deltaT_304_M106/scaling_factor);
scaledTF = TF.y/1.722592812;  %using magic number to match old TF to new
disp('WARNING! using magic number to match old TF to new')
scaledTF_mag = abs(scaledTF);
scaledTF_phase = angle(scaledTF);
scaledTF_phase = scaledTF_phase./pi.*180;
scaled_result_file = 'E:\c_test\Sasis Project\Python MRI TF Scripts\300_104_3T_Scaled_Interp.csv';
scaled_result_TF = [scaledTF_mag.',scaledTF_phase.'];
csvwrite(scaled_result_file,scaled_result_TF);
