function [] = myplot(lab,cad)
figure
subplot(2,1,1)
plot(lab.x,abs(lab.y),'b')
%legend('Matlab')
hold on
plot(cad.x,abs(cad.y1+1j*cad.y2),'m');
%legend('Semcad')
subplot(2,1,2)
plot(lab.x,angle(lab.y),'b')
hold on
plot(cad.x,angle(cad.y1+1j*cad.y2),'m');
