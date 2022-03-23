function dv=inflation(t,v)
    global den g ms mw vL kf t1;
    
    %密织物平面圆形伞plane circular
    D0=8; %伞衣名义直径
    C=0.8; %阻力系数
    lambda=1.74; %充气距离比例系数λ
    alpha=6; %Scheubel常数
    scale=0.04; %(CA)1/(CA)s

    %无收口non-close
    global tm;
    t1=lambda*D0/vL; %初始充气时间
    tm=alpha*D0/(0.9*vL); %充满时间(kv取0.9)
    kt=scale*C*(pi*(D0/2)^2)/t1;
    n=2;
    betat=(1-scale)*C*(pi*(D0/2)^2)/((tm-t1)^n); %n取2
      
    dv=zeros(4,1); 
    if t<=t1
        dv(1)=kt; %dCA/dt
    elseif t<=tm
        dv(1)=n*betat*(t-t1)^(n-1);
    else
        dv(1)=0;
    end
    dv(2)=1.5*den*kf*(v(1)^0.5)*dv(1); %dmf/dt
    dv(3)=(mw+ms)*g*sind(v(5))/(mw+ms+v(2))-0.5*den*(v(3)^2)*v(1)/(mw+ms+v(2))-v(3)*dv(2)/(mw+ms+v(2)); %dv/dt,缺(CA)w
    dv(4)=v(3); %ds/dt
    dv(5)=g*cosd(v(5))/v(3); %d(theta)/dt
    dv(6)=v(3)*cosd(v(5)); %dx/dt
    dv(7)=-v(3)*sind(v(5)); %dy/dt

end

