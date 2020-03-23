function ug = dufix(dh)
ug = 0;
for dh = -1:-1:dh
    ud = duda(1175+Re,dh)*deg;
    ug = ug + ud;
end