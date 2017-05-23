function barerrorbarExample

x = 0.2*randn(3,4,100) + 1;
xMeans = mean(x,3);
xMeansConf = repmat(2*0.2/10,3,4);
xMeansL = repmat(3*0.2/10,3,4);
xMeansU = repmat(4*0.2/10,3,4);

figure
barerrorbar(xMeans,xMeansConf);

figure
barerrorbar({3:5,xMeans,'m'}, {repmat((3:5)',1,4),xMeans, xMeansL,xMeansU,'bx'});

figure
barerrorbar({3:5,xMeans,'k'}, {repmat((3:5)',1,4),xMeans, xMeansL,xMeansU,'bx'});
hold on
barerrorbar({7:9,xMeans}, {repmat((7:9)',1,4),xMeans, 2*xMeansL,4*xMeansU,'d','color',[0.7 0.7 0.7]});
