clear
close all

epoch = [2020,9,1,4,0,0];
h0 = 1100;
a = h0 + Re;
kp = [a,0.001,85,120,60,0];
F = 0.015;
Isp = 1300;
dm = -F/Isp/9.8;
Mass = 200;
s(1,:) = mexOrbitDyn('sat','init',[epoch,kp,Mass]);
mexOrbitDyn('sat','SetForce',[8,1+2+4+8+16+32]);
mexOrbitDyn('sat','SetAD',5); 
mexOrbitDyn('sat','SetSRP',5);
m(1,:) = mexOrbitDyn('sat','GetMean');

%% ��ƹ��
kpstd = [Re + 1175,0.001,85,120,60,0];
sstd(1,:) = mexOrbitDyn('std','init',[epoch,kpstd,Mass]);
mexOrbitDyn('std','SetForce',[8,1+2+8+16]); % �����������ѹ
mstd(1,:) = mexOrbitDyn('std','GetMean');
%% 
step = 60;
i = 1;
%% �����������
ut = 15; % Ŀ����λ�ǣ���200~360��0~160��Χ��
ug = dufix(a - kpstd(1));
u1 = ut - ug - 0; % Ŀ����λ�Ǽ�ȥ�̶���λ�仯��
if u1 < 0
    u1 = u1 + 360;
end
du = duda(kpstd(1),a - kpstd(1))*deg; % ͣ�������λƯ����
dt = u1/du * a2t(a); % ͣ�����ͣ��ʱ��

%�ص��ʱ��,��һ��Ϊ��ʼʱ�䣬�ڶ���Ϊ����ʱ�䣬������������
% fireoff = [0,1440;...
%             1440*2,1440*3;...
%            1440*4,1440*5];
fireoff = [0,dt/step];
hpause = 1180; % �м�ͣ������Ĺ���߶�
for i = 2:(1440*10)
    qbi = c2q(getcoi(s(i-1,2:7)));
    if any(i>fireoff(:,1) & i<fireoff(:,2))        
        s(i,:) = mexOrbitDyn('sat','step',step);
    else
        s(i,:) = mexOrbitDyn('sat','step',[step, qbi', F, 0, 0, dm]);
    end
    sstd(i,:) = mexOrbitDyn('std','step',step);
    m(i,:) = mexOrbitDyn('sat','GetMean');
    mstd(i,:) = mexOrbitDyn('std','GetMean');
    if i > fireoff(end,2) && m(i,1) - Re > hpause
        hpause = 1180; % �´�ͣ������߶�
        u = ewmu(s(i,3),s(i,6),s(i,7));
        ustd = ewmu(sstd(i,3),sstd(i,6),sstd(i,7));
        ue = limitpi(u-ustd)*deg;
        ug = dufix(m(i,1) - kpstd(1));
        dltu = ut - ue - ug;
        if dltu < 0
            dltu = dltu + 360;
        end
        dotu = duda(kpstd(1),m(i,1) - kpstd(1))*deg;
        dt1 = dltu/dotu * a2t(m(i,1));
        fireoff = [fireoff;i,i+dt1/step];
    elseif m(i,1) - Re > 1175
        break;
    end
end
mexOrbitDyn('sat','remove');
mexOrbitDyn('std','remove');

t = (s(:,1)-s(1,1));
u = ewmu(s(:,3),s(:,6),s(:,7));
ustd = ewmu(sstd(:,3),sstd(:,6),sstd(:,7));
ue = limitpi(u-ustd)*deg;
ue(ue<0) = ue(ue<0)+360;
disp(ue(end));
disp(fireoff);
%%
figure('Name','˲����'),plot6coe(t,kdeg(s(:,2:7)));
figure('Name','ƽ����'),plot6coe(t,kdeg(m));
figure('Name','��λ��'),plot(t,ue);