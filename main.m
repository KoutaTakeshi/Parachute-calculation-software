global den g ms msh msy vZ hZ tm kf theta;
global C D0 A0 densh densy densystar Lsh Lxt b;
global kt betat n ti;
global t2 v2 Fk t2n v2n dtrailx dtraily itrailx itraily;
global sk FL me tmax;

me =0;
g = 9.8;
kf = 0.41; %附加质量系数kf = 0.752 / Ct^1.5
A0 = pi * (D0 / 2)^2; %伞衣面积
b = D0 / 100; %伞衣底边宽度
Lxt = Lsh + D0 / 2; %伞系统长度
msh = densh * Lsh; %伞绳质量
msy = densystar * b + densy * (D0 / 2 - b); %伞衣质量
ms = msh + msy; %降落伞质量

%% 拉直过程
tspan1 = (0 : 0.01 : 2); %求解时间范围
[t1, v1] = ode45('deployment', tspan1, [0 0 vZ vZ theta 0 hZ]);

%拉直力求解
tmax = 1;
while v1(tmax, 2) < Lxt
    tmax = tmax + 1;
end
FL = zeros([1, 1]);
for i = 1 : tmax - 1
    if v1(i, 2) < Lsh
        FL(i, 1) = densh * (v1(i, 3) - v1(i, 4))^2;
    elseif v1(i, 2) < Lxt
        FL(i, 1) = densy * (Lxt - v1(i, 2)) * (v1(i, 3) - v1(i, 4))^2; %伞衣单位长度质量
    end
end
dtrailx = [v1(1 : tmax - 1, 6) - v1(1 : tmax - 1, 2) .* cosd(v1(1 : tmax - 1, 5)), v1(1 : tmax - 1, 6)]; 
dtraily = [v1(1 : tmax - 1, 7) + v1(1 : tmax - 1, 2) .* sind(v1(1 : tmax - 1, 5)), v1(1 : tmax - 1, 7)];
v1(tmax, 4) = v1(tmax - 1, 3); %全长拉直时刻，引导伞速度突变为载荷速度

%% 充气过程
lambda = 1.74; %充气距离比例系数λ
lambda1 = 14;
lambdat = 4.42;
alpha = 6; %Scheubel常数
scale = 0.04; %(CA)1 / (CA)
vjsk = 30; %解除收口时物伞系统速度
n = 2;
if sk == 0
    %无收口
    ti = lambda * D0 / v1(tmax - 1, 3); %初始充气时间
    tm = alpha * D0 / (0.9 * v1(tmax - 1, 3)); %充满时间(kv取0.9)
    kt = scale * C * A0 / ti;
    betat = (1 - scale) * C * A0 / ((tm - ti)^n); %n取2
else
    %一次收口
    ti = lambda1 * scalesk * D0 / v1(tmax - 1, 3); %第一段充气时间
    tm = lambdat * (1 - scalesk) * (C * A0)^0.5 / vjsk; %充满时间
    kt = scalesk^2 * C * A0 / ti;
    betat = (C * A0 - kt * ti) / tm;
end

tspan2 = (0 : 0.01 : tm + 2.5);
[t2, v2] = ode45('inflation', tspan2, [0 0 v1(tmax - 1, 3) 0 v1(tmax - 1, 5) v1(tmax - 1, 6) v1(tmax - 1, 7)]);
[t2n, ~] = size(t2);
[v2n, ~] = size(v2);
itrailx = [v2(:, 6) - Lxt * cosd(v2(:, 5)), v2(:, 6)]; 
itraily = [v2(:, 7) + Lxt * sind(v2(:, 5)), v2(:, 7)];

%开伞动载求解
a1 = (ms + v2(1 : length(t2) - 1, 2)) .* diff(v2( : , 3)) ./ diff(t2);
a2 = 0.5 * den * (v2(1 : length(t2), 3).^2) .* v2(1 : length(t2), 1);
a3 = v2(1 : length(t2) - 1, 3) .* diff(v2( : , 2)) ./ diff(t2);
a1(length(t2)) =  a1(length(t2) - 1);
a3(length(t2)) = a3(length(t2) - 1);
Fk = a1 + a2 + a3 - ms * g * sind(v2(1 : length(t2), 5));

