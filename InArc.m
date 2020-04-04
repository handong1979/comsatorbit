% 判断幅角a是否在s~e的弧段内
function b = InArc(a,s,e)
a = mod(a,360);
s = mod(s,360);
e = mod(e,360);
if s < e
    if a > s && a < e
        b = 1;
    else
        b = 0;
    end
elseif s > e
    if a < e
        a = a + 360;
    end
    e = e + 360;
    if a > s && a < e
        b = 1;
    else
        b = 0;
    end
else
    b = 0;
end