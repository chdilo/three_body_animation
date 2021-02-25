function map = colorMap(color, m)
% 输出"黑→color色→白"渐变的颜色图
color = 0.5*color./(0.299*color(1)+0.587*color(2)+0.114*color(3));
cmax = max(color);
if cmax > 1
    color = color/cmax;
end
x = [0 0.5 1];
y1 = [0 color(1) 1];
y2 = [0 color(2) 1];
y3 = [0 color(3) 1];

p1 = polyfit(x, y1, 2);
p2 = polyfit(x, y2, 2);
p3 = polyfit(x, y3, 2);

map(:,1) = polyval(p1, linspace(0,1,m));
map(:,2) = polyval(p2, linspace(0,1,m));
map(:,3) = polyval(p3, linspace(0,1,m));

map(map>1) = 1;
map(map<0) = 0;
end
