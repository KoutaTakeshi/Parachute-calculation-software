function dv=chongqi(t,v)
    
    global den g ms mw vL kf;
    
    %密织物平面圆形伞plane circular
    D0=5; %伞衣名义直径
    C=0.8; %阻力系数
    kd=2; %开伞动载系数
    lambda=1.74; %充气距离比例系数λ
    alpha=6; %Scheubel常数
    scale=0.04; %(CA)1/(CA)s

    %无收口non-close
    t1=lambda*D0/vL; %初始充气时间
    tm=alpha*D0/vL; %充满时间（缺少速度修正系数kv）
    kt=scale*C*(pi*(D0/2)^2)/t1;
    n=2;
    betat=(1-scale)*C*(pi*(D0/2)^2)/((tm-t1)^n); %n取2
    %mf=den*kf*(ncCA)^1.5; %附加质量
      
    dv=zeros(2,1); %2或3
    if t<=t1
       dv(1)=kt; %dCA/dt
    else
       dv(1)=n*betat*(t-t1)^(n-1);
    end
    dv(2)=1.5*den*kf*(v(1)^0.5)*dv(1); %dmf/dt
    dv(3)=-(mw+ms)*g/(mw+ms+v(2))-0.5*den*(v(3)^2)*v(1)/(mw+ms+v(2))-v(3)*dv(2)/(mw+ms+v(2)); %dv/dt
    v(4)=(ms+v(2))*dv(3)+0.5*den*(v(3)^2)*v(1)+v(3)*dv(2)+ms*g; %Fk
    dv(5)=v(3); %ds/dt
%     if t<=t1
%         v(1)=kt*t;
%     else
%         v(1)=scale*C*(pi*(D0/2)^2)+betat*(t-t1)^2; %n取2
%     end
%     v(2)=den*kf*(v(1))^1.5; %mf
%     a=(mw+ms)/(mw+ms+v(2));
%     b=v(1)/(mw+ms+v(2)); %缺(CA)w
%     c=1/(mw+ms+v(2));
%     dv(3)=-a*g-0.5*den*v(2)*b-c*v(3)*dv(2); %dv/dt
%     dv(4)=v(3); %ds/dt
%     a=(mw+ms)/(mw+ms+den*kf*(CA(t))^1.5);
%     b=CA(t)/(mw+ms+den*kf*(CA(t))^1.5); %缺(CA)w
%     c=1/(mw+ms+den*kf*(CA(t))^1.5);
%     dv(1)=-a*g-0.5*den*CA(t)*b-c*v(1)*dCA(t); %dv/dt
%     dv(2)=v(1); %ds/dt

end

%伞阻力特征CA
% function a=CA(t)
% 
%     global t1 kt scale C D0 betat 
%     if t<=t1
%         a=kt*t;
%     else
%         a=scale*C*(pi*(D0/2)^2)+betat*(t-t1)^2; %n取2
%     end
%     
% end
