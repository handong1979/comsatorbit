% ����������
clear
clc
close all
% ��ʼʱ��
epoch = [2020,9,1,4,0,0];
% ��ʼ����߶�
h0 = 1175;
a = h0 + Re;
% ��ʼ�������
kp = [a,0.001,85,120,60,0];

F = 0.015; % ����
Isp = 1300; % �ȳ�
dm = -F/Isp/9.8; % ������
Mass = 200; % ����
% �����ʼ��
s(1,:) = mexOrbitDyn('sat','init',[epoch,kp,Mass]);
mexOrbitDyn('sat','SetForce',[8,1+2+4+8+16+32]); % 8*8�ף����´�����ѹ
mexOrbitDyn('sat','SetAD',5); 
mexOrbitDyn('sat','SetSRP',5);
m(1,:) = mexOrbitDyn('sat','GetMean'); % ƽ����

%% ƫ����ʸ��������
% ��ǰƫ����ʸ��
ex = m(1,2)*cos(m(1,5));
ey = m(1,2)*sin(m(1,5));
% Ŀ��ƫ����ʸ��
et = 0.001174; % Ŀ��ƫ����
wt = 90;    % Ŀ����ص����
ext = et*cosd(wt);
eyt = et*sind(wt);
% ƫ����ʸ����
dex = ex - ext;
dey = ey - eyt;
fc = atan2(dey,dex);
if fc < 0
    fc = fc + 2*pi;
end
% ������
df = 30*rad;
fdon = fc-df;
fdoff = fc+df;
fion = fc+pi-df;
fioff = fc+pi+df;
%% 
step = 60; % ���沽��
i = 1;
for i = 2:(1440*5)
    qbi = c2q(getcoi(s(i-1,2:7)));
    f = ma2ta(m(i-1,2),m(i-1,6));
    if InArc(f*deg,fdon*deg,fdoff*deg) % ���״̬
        s(i,:) = mexOrbitDyn('sat','step',[step, qbi', F, 0, 0, dm]);
    elseif InArc(f*deg,fion*deg,fioff*deg) % ���״̬
        s(i,:) = mexOrbitDyn('sat','step',[step, qbi', -F, 0, 0, dm]);
    else % �ǵ��״̬
        s(i,:) = mexOrbitDyn('sat','step',step);
    end
    m(i,:) = mexOrbitDyn('sat','GetMean'); % ƽ����
    ex = m(i,2)*cos(m(i,5));
    ey = m(i,2)*sin(m(i,5));
    dex = ex - ext;
    dey = ey - eyt;
    if sqrt(dex.^2+dey.^2) < 1e-4
        break;
    end
end
mexOrbitDyn('sat','remove');
%% ���ݱ���ͻ�ͼ
t = (s(:,1)-s(1,1)); %ʱ��
ex = m(:,2).*cos(m(:,5));
ey = m(:,2).*sin(m(:,5));
%% ����
figure('Name','˲����'),plot6coe(t,kdeg(s(:,2:7)));
saveas(gcf,'����������˲����.png')
figure('Name','ƽ����'),plot6coe(t,kdeg(m));
saveas(gcf,'����������ƽ����.png')
figure('Name','ƫ����ʸ��'),plot(ex,ey),grid on,hold on;
plot(ex(1),ey(1),'ro');plot(ex(end),ey(end),'rd');
plot(ext,eyt,'ko');
saveas(gcf,'ƫ����ʸ��.png');
