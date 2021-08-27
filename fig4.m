function fig4


fsiz = [0 0 .7 .4];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 2;
nc = 5;

h(3:5) = sim_conditioned_suppression(nr,nc,3:5,0);
h(8:10) = sim_partial_reinforcement(nr,nc,8:10,0);

for i = [1 6]
    h(i) = subplot(nr,nc,i);
    set(h(i),'visible','off');
end
h([2 7]) = [];

% --------
fn = def('fn');
fsA = def('fsA');
xsA = -.25;def('xsA');
ysA = def('ysA');
abc = def('abc');


for i= [2 3 4 6 7 8]
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
end


xsA = -1;
for i= [1 5]
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
end

end
