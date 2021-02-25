function dy = three_body(~, y, m)
% 三体的二阶微分方程
% m：质量

G = 6.67259e-11; % 万有引力常量
r1 = y(1:2);     % 位置矢量
r2 = y(3:4);
r3 = y(5:6);
dr1 = y(7:8);    % 位置矢量对时间的一阶导(速度)
dr2 = y(9:10);
dr3 = y(11:12);

% 加速度
d2r1 = -G*m(2)*(r1-r2)/norm(r1-r2)^3 - G*m(3)*(r1-r3)/norm(r1-r3)^3;
d2r2 = -G*m(3)*(r2-r3)/norm(r2-r3)^3 - G*m(1)*(r2-r1)/norm(r2-r1)^3;
d2r3 = -G*m(1)*(r3-r1)/norm(r3-r1)^3 - G*m(2)*(r3-r2)/norm(r3-r2)^3;

dy(1:2) = dr1;
dy(3:4) = dr2;
dy(5:6) = dr3;
dy(7:8) = d2r1;
dy(9:10) = d2r2;
dy(11:12) = d2r3;

dy = dy(:);
end