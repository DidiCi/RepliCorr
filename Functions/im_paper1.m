function im_paper1(name,width,height)
pos = get(gcf, 'papersize');

set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'papersize',[width,height])
set(gcf,'paperposition',[0,0,width,height])

print('-dtiff','-r1200',[name  '.tif'])
print('-dpng','-r1200',[name  '.png'])

end