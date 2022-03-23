clc;
clear;
global den g ms mw vL hL kf tm theta;
den=1.205; %空气密度
g=9.8;
ms=5; %降落伞质量
mw=70; %物体质量
vL=50; %拉直速度
hL=1500; %拉直阶段结束高度
kf=0.41;
theta=82; %倾斜角

%常微分方程求解(不考虑倾斜角)
tspan=(0:0.01:tm+2); %求解时间范围

[t,v]=ode45('inflation',tspan,[0 0 vL 0 theta 0 hL]);

hC=hL-v(1:length(t)-1,4); %充气阶段高度变化
x=(ms+v(1:length(t)-1,2)).*diff(v(:,3))./diff(t);
y=0.5*den*(v(1:length(t)-1,3).^2).*v(1:length(t)-1,1);
z=v(1:length(t)-1,3).*diff(v(:,2))./diff(t);
Fk=x+y+z-ms*g*sind(v(1:length(t)-1,5));

figure(1);

subplot(2,4,1);
plot(t,v(:,1)); 
title("CA");

subplot(2,4,2);
plot(t,v(:,2));
title("mf");

subplot(2,4,3);
plot(t,v(:,3));
title("v");

subplot(2,4,4);
plot(t,v(:,4));
title("s"); %位移

subplot(2,4,5);
plot(t,v(:,5));
title("theta");

subplot(2,4,6);
plot(t,v(:,6));
title("x"); %水平方向位移

subplot(2,4,7);
plot(t,v(:,7));
title("hC"); %竖直方向位移


figure(2);

plot(t(1:length(t)-1,:),Fk);
title("Fk");







