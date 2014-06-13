%説明用の図を書くスクリプト
streamline1=1.*[z',x'];
streamline2=2.*[z',x'];
streamline3=3.*[z',x'];

stt1 = [-1*cos(beta),-1*sin(beta)]
stt2 = [-2*cos(beta),-2*sin(beta)]
stt3 = [-3*cos(beta),-3*sin(beta)]

conex = [0.5;0;0.5];
coney = [0.5*tan(beta-delta);0;-0.5*tan(beta-delta)];
cone1 =[conex,coney]+stt1;
cone2 =[conex,coney]+stt2;
cone3 =[conex,coney]+stt3;

figure(4)
hold on
h1=plot(cone1(:,1),cone1(:,2))
set(h1,'linewidth',3)
h2=plot(cone2(:,1),cone2(:,2),'r')
set(h2,'linewidth',3)
h3=plot(cone3(:,1),cone3(:,2),'g')
set(h3,'linewidth',3)
h4=plot([0,stt3(1,1)],[0,stt3(1,2)],'k')
set(h4,'linewidth',3)
h5=plot(streamline1(:,1),streamline1(:,2))
set(h5,'linewidth',3)
h6=plot(streamline2(:,1),streamline2(:,2),'r')
set(h6,'linewidth',3)
h7=plot(streamline3(:,1),streamline3(:,2),'g')
set(h7,'linewidth',3)