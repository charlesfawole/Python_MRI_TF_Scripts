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
        result = abs((Etan.y*TF.y.')/100)
    case 'T'
        result = (abs((Etan.y*TF.y.')/100))^2
end
