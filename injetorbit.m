clear
clc
close all
% ��ʼʱ��
epoch = [2020,9,1,4,0,0];
% ��ʼ����߶�
h0 = 1100;
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

%% ��ƹ�����뵱ǰ�����ȣ�ֻ�а볤�᲻ͬ��
% ������λʱ���Ըñ�ƹ��Ϊ��׼
kpstd = [Re + 1170.28,0.001,85,120,60,0]; % TODO: ����ƽ������ת����˲����
% kpstd = [Re + 1170.622,0.001,85,120,60,0];
sstd(1,:) = mexOrbitDyn('std','init',[epoch,kpstd,Mass]);
mexOrbitDyn('std','SetForce',[8,1+2+8+16]); % 8*8�ף����£������������ѹ
mstd(1,:) = mexOrbitDyn('std','GetMean'); 
%% 
step = 60; % ���沽��
i = 1;
%% �����������
ut = 303; % Ŀ����λ��,����Զ��ڱ�ƹ������λ����203~360��0~160��Χ��(10��ת��)�� 
ug = dufix(a - kpstd(1)); % ��ǰ�����Ŀ�����������ʱ��׷�ϵ���λ��
u1 = ut - ug - 5; % Ŀ����λ�Ǽ�ȥ�̶���λ�仯��
if u1 < 0
    u1 = u1 + 360;
end
du = duda(kpstd(1),a - kpstd(1))*deg; % ��һͣ�������λƯ����
dt = u1/du * a2t(a); % ��һͣ�����ͣ��ʱ��

%�ص��ʱ��,��һ��Ϊ��ʼʱ��(step������)���ڶ���Ϊ����ʱ��(step������)��������������
% fireoff = [0,1440;...
%             1440*2,1440*3;...
%            1440*4,1440*5];
fireoff = [0,dt/step]; % ���ȳ�ʼ��Ϊ��һͣ�������ͣ��ʱ��
dh = 16; % ÿ����dh�߶ȣ�����һ��ͣ����������¼���ͣ��ʱ��
hpause = 1100 + dh; % �м�ͣ������Ĺ���߶�
for i = 2:(1440*10)
    qbi = c2q(getcoi(s(i-1,2:7)));
    if any(i>fireoff(:,1) & i<fireoff(:,2)) % ��������һ��ͣ�����
        s(i,:) = mexOrbitDyn('sat','step',step);
    else % ���������Ϊ���״̬
        s(i,:) = mexOrbitDyn('sat','step',[step, qbi', F, 0, 0, dm]);
    end
    sstd(i,:) = mexOrbitDyn('std','step',step); % ���Ƶõ�˲����
    m(i,:) = mexOrbitDyn('sat','GetMean'); % ƽ����
    mstd(i,:) = mexOrbitDyn('std','GetMean'); % ��ƹ��ƽ����
    if i > fireoff(end,2) && m(i,1) - Re > hpause %���ͣ�����ﵽ��һ��ͣ������߶�
        hpause = hpause + dh; % �´�ͣ������߶�
		% ������λƫ������������λ����ǰ��ȥ
        if hpause > 1175
            up = 0;
        else
            up = 5/75*(1175 - hpause); % ��λƫ�������߶Ȳ�����ȣ�75km��Ӧ5��
        end
        u = ewmu(m(i,2),m(i,5),m(i,6));
        ustd = ewmu(mstd(i,2),mstd(i,5),mstd(i,6));
        ue = limitpi(u-ustd)*deg; % ��ǰ��λ��
        ug = dufix(m(i,1) - mstd(1)); % �ӵ�ǰͣ���������ƹ���������ʱ��׷�ϵ���λ��
        dltu = ut - ue - ug - up; % ��ǰͣ�������Ҫ��������λ
        if dltu < 0
            dltu = dltu + 360;
        end
        if dltu > 360
            dltu = dltu - 360;
        end
        dotu = duda(mstd(i,1),m(i,1) - mstd(i,1))*deg; % ��ǰͣ�������λƯ����
        dt1 = dltu/dotu * a2t(m(i,1)); % ��ǰͣ�����ͣ��ʱ��
        fireoff = [fireoff;i,i+dt1/step]; % ��fireoff����������ӵ�ǰͣ�������ͣ��ʱ��
    elseif m(i,1) - Re > 1175 % �ﵽ��ƹ����ֹͣ
        break;
    end
end
mexOrbitDyn('sat','remove');
mexOrbitDyn('std','remove');

%% ���ݱ���ͻ�ͼ
t = (s(:,1)-s(1,1)); %ʱ��
u = ewmu(s(:,3),s(:,6),s(:,7));
ustd = ewmu(sstd(:,3),sstd(:,6),sstd(:,7));
ue = limitpi(u-ustd)*deg;
ue(ue<0) = ue(ue<0)+360; % ��λ��0~360

dhl = (1100:dh:1175)'; % ����ͣ������ĸ߶�
diary('diary.txt');
fprintf('Ŀ����λ��%f, ת�ƽ�����λ��%f\n',ut,ue(end));
fprintf('ͣ������߶�%.0fkm, ͣ��ʱ�� %.1f�� -- %.1f��,ͣ��ʱ�� %.1f��\n',[dhl,fireoff,fireoff(:,2)-fireoff(:,1)]');
fprintf('ת��ʱ��%.2f��\n',i*step/86400);
diary off
s(:,2) = s(:,2) - Re; % �볤�ᶼ��ȥ����뾶
m(:,1) = m(:,1) - Re; % �볤�ᶼ��ȥ����뾶
%% ����
figure('Name','˲����'),plot6coe(t,kdeg(s(:,2:7)));
saveas(gcf,'˲����.png')
figure('Name','ƽ����'),plot6coe(t,kdeg(m));
saveas(gcf,'ƽ����.png')
figure('Name','��λ��'),plot(t,ue);
xlabel('t(��)'),ylabel('��λ(��)');
saveas(gcf,'��λ��.png')