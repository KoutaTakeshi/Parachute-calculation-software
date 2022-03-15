global den g ms mw vL kf;
den=1.205; %空气密度
g=9.8;
ms=1; %降落伞质量
mw=1; %物体质量
vL=50; %拉直速度
kf=0.41;
 
%常微分方程求解(不考虑倾斜角)
tspan=[0 100];

[t,v]=ode45('chongqi',tspan,[0 0 vL 0 0]);

plot(t,v(:,4)); %v图像
%plot(t,s);




