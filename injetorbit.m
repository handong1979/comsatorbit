clear
clc
close all
% 初始时间
epoch = [2020,9,1,4,0,0];
% 初始轨道高度
h0 = 1100;
a = h0 + Re;
% 初始轨道根数
kp = [a,0.001,85,120,60,0];

F = 0.015; % 推力
Isp = 1300; % 比冲
dm = -F/Isp/9.8; % 秒流量
Mass = 200; % 质量
% 轨道初始化
s(1,:) = mexOrbitDyn('sat','init',[epoch,kp,Mass]);
mexOrbitDyn('sat','SetForce',[8,1+2+4+8+16+32]); % 8*8阶，日月大气光压
mexOrbitDyn('sat','SetAD',5); 
mexOrbitDyn('sat','SetSRP',5);
m(1,:) = mexOrbitDyn('sat','GetMean'); % 平根数

%% 标称轨道，与当前轨道相比，只有半长轴不同。
% 计算相位时，以该标称轨道为基准
kpstd = [Re + 1170.28,0.001,85,120,60,0]; % TODO: 设置平根数，转化成瞬根数
% kpstd = [Re + 1170.622,0.001,85,120,60,0];
sstd(1,:) = mexOrbitDyn('std','init',[epoch,kpstd,Mass]);
mexOrbitDyn('std','SetForce',[8,1+2+8+16]); % 8*8阶，日月，不含大气与光压
mstd(1,:) = mexOrbitDyn('std','GetMean'); 
%% 
step = 60; % 仿真步长
i = 1;
%% 调相参数计算
ut = 303; % 目标相位角,是相对对于标称轨道的相位，在203~360或0~160范围内(10天转移)。 
ug = dufix(a - kpstd(1)); % 当前轨道到目标轨道连续变轨时，追赶的相位差
u1 = ut - ug - 5; % 目标相位角减去固定相位变化量
if u1 < 0
    u1 = u1 + 360;
end
du = duda(kpstd(1),a - kpstd(1))*deg; % 第一停泊轨道相位漂移率
dt = u1/du * a2t(a); % 第一停泊轨道停泊时间

%关点火时段,第一列为开始时间(step，分钟)，第二列为结束时间(step，分钟)，可以任意行数
% fireoff = [0,1440;...
%             1440*2,1440*3;...
%            1440*4,1440*5];
fireoff = [0,dt/step]; % 首先初始化为第一停泊轨道的停泊时间
dh = 16; % 每升高dh高度，增加一个停泊轨道，重新计算停泊时间
hpause = 1100 + dh; % 中间停泊轨道的轨道高度
for i = 2:(1440*10)
    qbi = c2q(getcoi(s(i-1,2:7)));
    if any(i>fireoff(:,1) & i<fireoff(:,2)) % 满足任意一个停泊轨道
        s(i,:) = mexOrbitDyn('sat','step',step);
    else % 其他情况均为点火状态
        s(i,:) = mexOrbitDyn('sat','step',[step, qbi', F, 0, 0, dm]);
    end
    sstd(i,:) = mexOrbitDyn('std','step',step); % 外推得到瞬根数
    m(i,:) = mexOrbitDyn('sat','GetMean'); % 平根数
    mstd(i,:) = mexOrbitDyn('std','GetMean'); % 标称轨道平根数
    if i > fireoff(end,2) && m(i,1) - Re > hpause %完成停泊、达到下一个停泊轨道高度
        hpause = hpause + dh; % 下次停泊轨道高度
		% 设置相位偏置量，避免相位超到前面去
        if hpause > 1175
            up = 0;
        else
            up = 5/75*(1175 - hpause); % 相位偏置与轨道高度差成正比，75km对应5度
        end
        u = ewmu(m(i,2),m(i,5),m(i,6));
        ustd = ewmu(mstd(i,2),mstd(i,5),mstd(i,6));
        ue = limitpi(u-ustd)*deg; % 当前相位差
        ug = dufix(m(i,1) - mstd(1)); % 从当前停泊轨道到标称轨道连续点火时，追赶的相位差
        dltu = ut - ue - ug - up; % 当前停泊轨道需要调整的相位
        if dltu < 0
            dltu = dltu + 360;
        end
        if dltu > 360
            dltu = dltu - 360;
        end
        dotu = duda(mstd(i,1),m(i,1) - mstd(i,1))*deg; % 当前停泊轨道相位漂移率
        dt1 = dltu/dotu * a2t(m(i,1)); % 当前停泊轨道停泊时间
        fireoff = [fireoff;i,i+dt1/step]; % 在fireoff数组最后增加当前停泊轨道的停泊时间
    elseif m(i,1) - Re > 1175 % 达到标称轨道后停止
        break;
    end
end
mexOrbitDyn('sat','remove');
mexOrbitDyn('std','remove');

%% 数据保存和画图
t = (s(:,1)-s(1,1)); %时间
u = ewmu(s(:,3),s(:,6),s(:,7));
ustd = ewmu(sstd(:,3),sstd(:,6),sstd(:,7));
ue = limitpi(u-ustd)*deg;
ue(ue<0) = ue(ue<0)+360; % 相位差0~360

dhl = (1100:dh:1175)'; % 各个停泊轨道的高度
diary('diary.txt');
fprintf('目标相位角%f, 转移结束相位角%f\n',ut,ue(end));
fprintf('停泊轨道高度%.0fkm, 停泊时间 %.1f分 -- %.1f分,停泊时长 %.1f分\n',[dhl,fireoff,fireoff(:,2)-fireoff(:,1)]');
fprintf('转移时间%.2f天\n',i*step/86400);
diary off
s(:,2) = s(:,2) - Re; % 半长轴都减去地球半径
m(:,1) = m(:,1) - Re; % 半长轴都减去地球半径
%% 曲线
figure('Name','瞬根数'),plot6coe(t,kdeg(s(:,2:7)));
saveas(gcf,'瞬根数.png')
figure('Name','平根数'),plot6coe(t,kdeg(m));
saveas(gcf,'平根数.png')
figure('Name','相位差'),plot(t,ue);
xlabel('t(天)'),ylabel('相位(度)');
saveas(gcf,'相位差.png')