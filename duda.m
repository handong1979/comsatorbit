% �볤���Сƫ�����ÿ��������ڵ���λƯ����(rad)
function du = duda(a,da)
n = sqrt(GEarth/a^3);
dn = -1.5*da/a*n;
du = dn*2*pi/n;