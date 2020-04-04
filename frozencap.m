% 冻结轨道捕获
clear
clc
close all
% 初始时间
epoch = [2020,9,1,4,0,0];
% 初始轨道高度
h0 = 1175;
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

%% 偏心率矢量调整量
% 当前偏心率矢量
ex = m(1,2)*cos(m(1,5));
ey = m(1,2)*sin(m(1,5));
% 目标偏心率矢量
et = 0.001174; % 目标偏心率
wt = 90;    % 目标近地点幅角
ext = et*cosd(wt);
eyt = et*sind(wt);
% 偏心率矢量差
dex = ex - ext;
dey = ey - eyt;
fc = atan2(dey,dex);
if fc < 0
    fc = fc + 2*pi;
end
% 点火参数
df = 30*rad;
fdon = fc-df;
fdoff = fc+df;
fion = fc+pi-df;
fioff = fc+pi+df;
%% 
step = 60; % 仿真步长
i = 1;
for i = 2:(1440*5)
    qbi = c2q(getcoi(s(i-1,2:7)));
    f = ma2ta(m(i-1,2),m(i-1,6));
    if InArc(f*deg,fdon*deg,fdoff*deg) % 点火状态
        s(i,:) = mexOrbitDyn('sat','step',[step, qbi', F, 0, 0, dm]);
    elseif InArc(f*deg,fion*deg,fioff*deg) % 点火状态
        s(i,:) = mexOrbitDyn('sat','step',[step, qbi', -F, 0, 0, dm]);
    else % 非点火状态
        s(i,:) = mexOrbitDyn('sat','step',step);
    end
    m(i,:) = mexOrbitDyn('sat','GetMean'); % 平根数
    ex = m(i,2)*cos(m(i,5));
    ey = m(i,2)*sin(m(i,5));
    dex = ex - ext;
    dey = ey - eyt;
    if sqrt(dex.^2+dey.^2) < 1e-4
        break;
    end
end
mexOrbitDyn('sat','remove');
%% 数据保存和画图
t = (s(:,1)-s(1,1)); %时间
ex = m(:,2).*cos(m(:,5));
ey = m(:,2).*sin(m(:,5));
%% 曲线
figure('Name','瞬根数'),plot6coe(t,kdeg(s(:,2:7)));
saveas(gcf,'冻结轨道捕获瞬根数.png')
figure('Name','平根数'),plot6coe(t,kdeg(m));
saveas(gcf,'冻结轨道捕获平根数.png')
figure('Name','偏心率矢量'),plot(ex,ey),grid on,hold on;
plot(ex(1),ey(1),'ro');plot(ex(end),ey(end),'rd');
plot(ext,eyt,'ko');
saveas(gcf,'偏心率矢量.png');
