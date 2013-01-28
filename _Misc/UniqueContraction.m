clc;
clear('all');

x = randi(10,1,20);
[xuni,m,n] = unique(x);

disp(x.^2);
disp(xuni(n).^2);
