clear, clc, close all
%% 求解三体的微分方程
MO = 2e30;              % 太阳质量
pc = 3e16;              % 秒差距
yr = 31557600;          % 儒略年
m = [0.75 1 1.25] * MO; % 质量
K = 6000 * sqrt(m/MO);  % 温度
L = (m/MO).^3.5;        % 光度
r0 = [1 3; -2 -1; 1 -1] * pc; % 初始位置
v0 = [0 0; 0 0; 0 0];      % 初速度
frameRate = 30;            % 帧率，Hz，FPS
n = 16;                    % 每帧用n张图平均
steps = n * frameRate;     % 每秒步数
dur = 2*60;                % 持续时间，second
nFrames = dur * frameRate; % 总帧数
xSpeed = 2e7 * yr;         % 倍速
time = (0:1/steps:dur-1/steps) * xSpeed;

% 指定绝对误差容限和相对误差容限
options = odeset('AbsTol', 1e-50, 'RelTol', 1e-13);
% 求解三体的二阶微分方程
[T, Y] = ode15s(@(t, x) three_body(t, x, m), ...
    time, [r0(1,:) r0(2,:) r0(3,:) v0(1,:) v0(2,:) v0(3,:)], options);

rx = Y(:, [1 3 5]);
ry = Y(:, [2 4 6]);
% vx = Y(:, [7 9 11]);
% vy = Y(:, [8 10 12]);
% G = 6.67259e-11;
% a(:,1:2) = -G*m(2)*(Y(:,1:2)-Y(:,3:4))./(sum(abs(Y(:,1:2)-Y(:,3:4)).^2,2).^(1/2)).^3 - G*m(3)*(Y(:,1:2)-Y(:,5:6))./(sum(abs(Y(:,1:2)-Y(:,5:6)).^2,2).^(1/2)).^3;
% a(:,3:4) = -G*m(3)*(Y(:,3:4)-Y(:,5:6))./(sum(abs(Y(:,3:4)-Y(:,5:6)).^2,2).^(1/2)).^3 - G*m(1)*(Y(:,3:4)-Y(:,1:2))./(sum(abs(Y(:,3:4)-Y(:,1:2)).^2,2).^(1/2)).^3;
% a(:,5:6) = -G*m(1)*(Y(:,5:6)-Y(:,1:2))./(sum(abs(Y(:,5:6)-Y(:,1:2)).^2,2).^(1/2)).^3 - G*m(2)*(Y(:,5:6)-Y(:,3:4))./(sum(abs(Y(:,5:6)-Y(:,3:4)).^2,2).^(1/2)).^3;
% ax = a(:,[1 3 5]);
% ay = a(:,[2 4 6]);
% 
% subplot(131)
% plot(rx(:,1),ry(:,1),rx(:,2),ry(:,2),rx(:,3),ry(:,3))
% title('位置矢量')
% grid on
% axis image
% axis([-4*16/9 4*16/9 -4 4]*pc)
% subplot(132)
% plot(vx(:,1),vy(:,1),vx(:,2),vy(:,2),vx(:,3),vy(:,3))
% title('速度矢量')
% grid on
% axis image
% subplot(133)
% plot(ax(:,1),ay(:,1),ax(:,2),ay(:,2),ax(:,3),ay(:,3))
% title('加速度矢量')
% grid on
% axis image
% 
% ra(:,1) = sqrt(Y(:,1).^2 + Y(:,2).^2);
% ra(:,2) = sqrt(Y(:,3).^2 + Y(:,4).^2);
% ra(:,3) = sqrt(Y(:,5).^2 + Y(:,6).^2);
% va(:,1) = sqrt(Y(:,7).^2 + Y(:,8).^2);
% va(:,2) = sqrt(Y(:,9).^2 + Y(:,10).^2);
% va(:,3) = sqrt(Y(:,11).^2 + Y(:,12).^2);
% aa(:,1) = sqrt(a(:,1).^2 + a(:,2).^2);
% aa(:,2) = sqrt(a(:,3).^2 + a(:,4).^2);
% aa(:,3) = sqrt(a(:,5).^2 + a(:,6).^2);
% figure
% subplot(311)
% plot(T, ra)
% title('位置矢量的模')
% subplot(312)
% plot(T, va)
% title('速度矢量的模')
% subplot(313)
% plot(T, aa)
% title('加速度矢量的模')

%% 生成动画帧序列
width = 640;
height = 360;
yLim = [-4 4]*pc;
xLim = width/height*[-4 4]*pc;

colorMap1 = colorMap(colorTemp10deg(K(1)), 65536);
colorMap2 = colorMap(colorTemp10deg(K(2)), 65536);
colorMap3 = colorMap(colorTemp10deg(K(3)), 65536);

radius = height/256;
brightness = (height/radius/100)^2 * L;
ringRadius = height/8;

rx = (rx-xLim(1))/(xLim(2)-xLim(1)) * width;
ry = (ry-yLim(1))/(yLim(2)-yLim(1)) * height;
ry = height - ry;

xSequence = gpuArray(1:width);
ySequence = gpuArray(1:height);

path = uigetdir('.\', '选择动画帧序列保存路径');
Fig = waitbar(0, '正在处理帧...');
for i = 1:nFrames
    image1 = zeros(height, width, 'gpuArray');
    image2 = zeros(height, width, 'gpuArray');
    image3 = zeros(height, width, 'gpuArray');
    for j = 1:n
        [xArray1, yArray1] = meshgrid(xSequence-rx((i-1)*n+j,1), ySequence-ry((i-1)*n+j,1));
        [xArray2, yArray2] = meshgrid(xSequence-rx((i-1)*n+j,2), ySequence-ry((i-1)*n+j,2));
        [xArray3, yArray3] = meshgrid(xSequence-rx((i-1)*n+j,3), ySequence-ry((i-1)*n+j,3));
        imageTemp1 = brightness(1)./(1 + (xArray1/radius).^2 + (yArray1/radius).^2) + ...
            brightness(1)/2./(1 + (xArray1/radius/8).^2 + (4*yArray1/radius).^2) + ...
            brightness(1)/64./(1 + ((sqrt(xArray1.^2 + yArray1.^2)-ringRadius)/(4*radius)).^2) + ...
            brightness(1)/64./(1 + ((sqrt(xArray1.^2 + yArray1.^2)+ringRadius)/(4*radius)).^2);
        imageTemp2 = brightness(2)./(1 + (xArray2/radius).^2 + (yArray2/radius).^2) + ...
            brightness(2)/2./(1 + (xArray2/radius/8).^2 + (4*yArray2/radius).^2) + ...
            brightness(2)/64./(1 + ((sqrt(xArray2.^2 + yArray2.^2)-ringRadius)/(4*radius)).^2) + ...
            brightness(2)/64./(1 + ((sqrt(xArray2.^2 + yArray2.^2)+ringRadius)/(4*radius)).^2);
        imageTemp3 = brightness(3)./(1 + (xArray3/radius).^2 + (yArray3/radius).^2) + ...
            brightness(3)/2./(1 + (xArray3/radius/8).^2 + (4*yArray3/radius).^2) + ...
            brightness(3)/64./(1 + ((sqrt(xArray3.^2 + yArray3.^2)-ringRadius)/(4*radius)).^2) + ...
            brightness(3)/64./(1 + ((sqrt(xArray3.^2 + yArray3.^2)+ringRadius)/(4*radius)).^2);
        image1 = image1 + imageTemp1/n;
        image2 = image2 + imageTemp2/n;
        image3 = image3 + imageTemp3/n;
    end
    rgbImage = ind2rgb(uint16(65536*image1), colorMap1) + ...
        ind2rgb(uint16(65536*image2), colorMap2)+...
        ind2rgb(uint16(65536*image3), colorMap3);
    imwrite(rgbImage, fullfile(path, sprintf('%06u.png', i)))
    waitbar(i/nFrames, Fig, ...
        sprintf('正在处理帧...%.2f%% (%u/%u)', 100*i/nFrames, i, nFrames));
end
close(Fig)

%% 生成轨迹帧序列
image1 = zeros(height, width, 'gpuArray');
image2 = zeros(height, width, 'gpuArray');
image3 = zeros(height, width, 'gpuArray');

radius = height/1024;
brightness = (height/radius/(512*n))^2 * L;
halfLife = 8;
decayRate = 0.5^(1/(halfLife*steps));

path = uigetdir('.\', '选择轨迹帧序列保存路径');
Fig = waitbar(0, '正在处理帧...');
for i = 1:nFrames
    for j = 1:n
        [xArray1, yArray1] = meshgrid(xSequence-rx((i-1)*n+j,1), ySequence-ry((i-1)*n+j,1));
        [xArray2, yArray2] = meshgrid(xSequence-rx((i-1)*n+j,2), ySequence-ry((i-1)*n+j,2));
        [xArray3, yArray3] = meshgrid(xSequence-rx((i-1)*n+j,3), ySequence-ry((i-1)*n+j,3));
        image1 = decayRate*image1 + brightness(1)./(1 + (xArray1/radius).^2 + (yArray1/radius).^2);
        image2 = decayRate*image2 + brightness(2)./(1 + (xArray2/radius).^2 + (yArray2/radius).^2);
        image3 = decayRate*image3 + brightness(3)./(1 + (xArray3/radius).^2 + (yArray3/radius).^2);
    end
    rgbImage = ind2rgb(uint16(65536*image1), colorMap1) + ...
        ind2rgb(uint16(65536*image2), colorMap2) + ...
        ind2rgb(uint16(65536*image3), colorMap3);
    imwrite(rgbImage, fullfile(path, sprintf('%06u.png', i)))
    waitbar(i/nFrames, Fig, ...
        sprintf('正在处理帧...%.2f%% (%u/%u)', 100*i/nFrames, i, nFrames));
end
close(Fig)

%% 生成音频
[file, path] = uiputfile({'*.flac';'*.wav'}, '保存音频文件', 'Audio');
fs = 48e3; % 音频采样率
t = 0:1/fs:dur-1/fs;
kf = 16384; % 调频灵敏度
fc = 128;   % 载波频率

x1 = Y(:,1); % x轴坐标
x2 = Y(:,3);
x3 = Y(:,5);

v1 = sqrt(Y(:,7).^2 + Y(:,8).^2); % 速度
v2 = sqrt(Y(:,9).^2 + Y(:,10).^2);
v3 = sqrt(Y(:,11).^2 + Y(:,12).^2);

T = T / xSpeed;
f1 = interp1(T, v1, t, 'spline');
f2 = interp1(T, v2, t, 'spline');
f3 = interp1(T, v3, t, 'spline');

fMax = max([f1 f2 f3]);
f1 = f1 / fMax;
f2 = f2 / fMax;
f3 = f3 / fMax;

fIntegral1 = cumtrapz(t, f1);
fIntegral2 = cumtrapz(t, f2);
fIntegral3 = cumtrapz(t, f3);

lA(:,1) = (xLim(2)-x1) / (xLim(2)-xLim(1));
lA(:,2) = (xLim(2)-x2) / (xLim(2)-xLim(1));
lA(:,3) = (xLim(2)-x3) / (xLim(2)-xLim(1));
lA(lA<0) = 0;
lA(lA>1) = 1;
leftAmplitude1 = interp1(T, lA(:,1), t, 'spline');
leftAudio1 = f1.*leftAmplitude1.*cos(2*pi * (fc*t + kf*fIntegral1));
leftAmplitude2 = interp1(T, lA(:,2), t, 'spline');
leftAudio2 = f2.*leftAmplitude2.*cos(2*pi * (fc*t + kf*fIntegral2));
leftAmplitude3 = interp1(T, lA(:,3), t, 'spline');
leftAudio3 = f3.*leftAmplitude3.*cos(2*pi * (fc*t + kf*fIntegral3));

leftAudio = leftAudio1 + leftAudio2 + leftAudio3;

rA(:,1) = (x1-xLim(1)) / (xLim(2)-xLim(1));
rA(:,2) = (x2-xLim(1)) / (xLim(2)-xLim(1));
rA(:,3) = (x3-xLim(1)) / (xLim(2)-xLim(1));
rA(rA<0) = 0;
rA(rA>1) = 1;
rightAmplitude1 = interp1(T, rA(:,1), t, 'spline');
rightAudio1 = f1.*rightAmplitude1.*cos(2*pi * (fc*t + kf*fIntegral1));
rightAmplitude2 = interp1(T, rA(:,2), t, 'spline');
rightAudio2 = f2.*rightAmplitude2.*cos(2*pi * (fc*t + kf*fIntegral2));
rightAmplitude3 = interp1(T, rA(:,3), t, 'spline');
rightAudio3 = f3.*rightAmplitude3.*cos(2*pi * (fc*t + kf*fIntegral3));

rightAudio = rightAudio1 + rightAudio2 + rightAudio3;

y = [leftAudio' rightAudio'];
y = y/max(abs(y), [], 'all');
audiowrite([path, file], y, fs, 'BitsPerSample', 24)
