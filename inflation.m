function dv = inflation(t, v)
    global den g ms mw kf tm;
    global kt betat n ti sk ts CAw;
    
    dv = zeros(7, 1);
    if sk == 0 %无收口
        if t <= ti
            dv(1) = kt; %dCA / dt
        elseif t <= tm
            dv(1) = n * betat * (t - ti).^(n - 1);
        else
            dv(1) = 0;
        end
    else %一次收口
        if t <= ti
            dv(1) = kt; %dCA / dt
        elseif t <= (ti + ts)
            dv(1) = 0;
        elseif t <= (ti + ts + tm)
            dv(1) = 2 * betat * (t - ti - ts);
        else
            dv(1) = 0;
        end
    end
    dv(2) = 1.5 * den * kf * (v(1)^0.5) * dv(1); %dmf / dt
    dv(3) = (mw + ms) * g * sind(v(5)) / (mw + ms + v(2)) - 0.5 * den * (v(3)^2) * (v(1) + CAw) / (mw + ms + v(2))...
                 - v(3) * dv(2) / (mw + ms + v(2)); %dv / dt
    dv(4) = v(3); %ds / dt
    dv(5) = (180 / pi) * g * cosd(v(5)) / v(3); %d(theta) / dt
    dv(6) = v(3) * cosd(v(5)); %dx / dt
    dv(7) =  - v(3) * sind(v(5)); %dh / dt

end
