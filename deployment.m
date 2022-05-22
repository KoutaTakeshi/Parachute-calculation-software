function dv = deployment(t, v)
    global g mw den msh msy densh densy densystar D0 Lsh Lxt CAys mys b CAw;
    
    dv = zeros(7, 1);
    
    dv(2) = v(3) - v(4); %dL / dt
    if v(2) <= Lsh
        dv(1) = densh * dv(2); %dme / dt(me:已拉出部分质量)
        dv(3) = g * sind(v(5)) - ((0.5 * den * v(3)^2 * CAw) + densh * ((v(3) - v(4))^2)) / (mw + v(1)); %dvw / dt(vw:物体速度)
    elseif v(2) <= (Lsh+b)
        dv(1) = densystar * dv(2);
        dv(3) = g * sind(v(5)) - ((0.5 * den * v(3)^2 * CAw) + densystar * (Lxt - v(2)) * ((v(3) - v(4))^2)) / (mw + v(1));
    elseif v(2) <= Lxt
        dv(1) = (msy - densystar * b) * 8 * (Lxt - v(2)) / (D0^2) * dv(2);
        dv(3) = g * sind(v(5)) - ((0.5 * den * v(3)^2 * CAw) + densy * (Lxt - v(2)) * ((v(3) - v(4))^2)) / (mw + v(1));
    end 
    dv(4) = g * sind(v(5)) - (0.5 * den * (v(4)^2) * CAys) / (mys + (msh + msy - v(1))); 
    %dvys / dt(vys:引导伞速度),忽略Fd(伞袋阻力), Fsh、Fsy(伞绳、伞衣拉出阻力)
    dv(5) = (180 / pi) * g * cosd(v(5)) / v(3); %d(theta) / dt
    dv(6) = v(3) * cosd(v(5)); %dx / dt
    dv(7) =  - v(3) * sind(v(5)); %dy / dt
    
end
    
    