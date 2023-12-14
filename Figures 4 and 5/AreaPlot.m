function h=AreaPlot(t,y,SE,LineColor,Transparency,LWide)
% This is used for ploting the transparency figures.
% Input: t: time
%        y: mean value
%        SE: standord error
% Author: Xiong Xiao; Date: 2013-07-12

%%
h=plot(t,y,'Color',LineColor,'LineWidth',LWide);
hold on
upper=y+SE;
lower=y-SE;
yaxis=[upper,fliplr(lower)];
xaxis=[t,fliplr(t)];
fill(xaxis,yaxis,LineColor,'edgecolor','None');
alpha(Transparency);

end