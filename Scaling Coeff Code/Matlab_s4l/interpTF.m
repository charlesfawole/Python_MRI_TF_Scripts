TF.position = flipud(xlsread('104 to 106.xlsx','1.5T','A2:A37'));
TF.mag0 = xlsread('104 to 106.xlsx','1.5T','E2:E37');
TF.phase0 = xlsread('104 to 106.xlsx','1.5T','I2:I37');
TF.complex = TF.mag0.*exp(j*((TF.phase0/180)*pi));

TF.y = interp1(TF.position,TF.complex,(1:1:43),'linear','extrap');
TF.x = 0.01:0.01:0.43;
TF.mag = abs(TF.y);
TF.phase = angle(TF.y);
TF.phase = TF.phase./pi.*180;
result_file = 'interpTF104_1_5T.csv';
TF_result = [TF.mag.',TF.phase.'];
csvwrite(result_file,TF_result);
%xlswrite(result_file,TF_result,strcat('A1',':','B',num2str(length(TF.phase))));

figure
subplot(2,1,1)
plot(TF.x*100,TF.mag,'b')
hold on
plot(TF.position,TF.mag0,'r')
hold off
title('scaled transfer function magnitude')
subplot(2,1,2)
plot(TF.x*100,TF.phase,'b');
hold on
plot(TF.position,TF.phase0,'r')
hold off
title('scaled transfer function phase')
