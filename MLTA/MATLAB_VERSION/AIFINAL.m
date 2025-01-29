clear all
close all
clc
warning off

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 12)


% Dataset
A = load('DS.txt'); % (fi, H [A/m], omega [Hz], a [m]) and  (Tc(°C), t(s))
A(:,1) = A(:,1)/10; % percentagem fixed
X = A(:,1:4);
Y = A(:,5:6);

c = cvpartition(size(Y,1),"Holdout",0.2);
XTrain = X(training(c),:);
XTest = X(test(c),:);
YTrain1 = Y(training(c),1);
YTest1 = Y(test(c),1);
YTrain2 = Y(training(c),2);
YTest2 = Y(test(c),2);

    Mdl1 = fitrnet(...
    XTrain, ...
    YTrain1, ...
    'LayerSizes', [161 10 300], ...
    'Activations', 'sigmoid', ...
    'Lambda', 5.0754e-06, ...
    'IterationLimit', 1000, ...
    'Standardize', true);


    Mdl2 = fitrnet(...
    XTrain, ...
    YTrain2, ...
    'LayerSizes', [10 293], ...
    'Activations', 'sigmoid', ...
    'Lambda', 2.1568e-08, ...
    'IterationLimit', 1000, ...
    'Standardize', true);


testPredictions1 = predict(Mdl1,XTest);
testPredictions2 = predict(Mdl2,XTest);

figure,
plot(YTest1,testPredictions1,".")
hold on
plot(YTest1,YTest1)
hold off
xlabel("True Temperature (°C)")
ylabel("Predicted Temperature (°C)")
print(gcf,'Temperature.png','-dpng','-r300');
R2adj1 = R2A(YTest1,testPredictions1,size(X,2))

figure,
plot(YTest2,testPredictions2,".")
hold on
plot(YTest2,YTest2)
hold off
xlabel("True Time (s)")
ylabel("Predicted Time (s)")
print(gcf,'Time.png','-dpng','-r300');
R2adj2 = R2A(YTest2,testPredictions2,size(X,2))

R1 = [XTest YTest1 testPredictions1];
R2 = [XTest YTest2 testPredictions2];

% Temperature Surfaces
x1lin = linspace(min(R1(:,1)),max(R1(:,1)),50);
y1lin = linspace(min(R1(:,2)),max(R1(:,2)),50);
[X1,Y1] = meshgrid(x1lin, y1lin);
Z1 = griddata(R1(:,1),R1(:,2),R1(:,end),X1,Y1,'v4'); %predict
Z2 = griddata(R1(:,1),R1(:,2),R1(:,5),X1,Y1,'v4'); %real
figure
surf(X1,Y1,Z1,'FaceAlpha',0.3,'FaceColor','r','LineWidth',0.1,'EdgeColor','r')
hold on
surf(X1,Y1,Z2,'FaceAlpha',0.3,'FaceColor','b','LineWidth',0.1,'EdgeColor','b')
%legend('ANN', 'Real');
xlabel('$\phi$ (\%)','interpreter','latex');
ylabel('$H_0$ (A/m)','interpreter','latex');
zlabel('$T_c$ ($^{\circ}$C)','interpreter','latex');
view([87.5 10])
grid off
box on
print(gcf,'1.png','-dpng','-r300');

x1lin = linspace(min(R1(:,3)),max(R1(:,3)),50);
y1lin = linspace(min(R1(:,4)),max(R1(:,4)),50);
[X1,Y1] = meshgrid(x1lin, y1lin);
Z1 = griddata(R1(:,3),R1(:,4),R1(:,end),X1,Y1,'v4'); %predict
Z2 = griddata(R1(:,3),R1(:,4),R1(:,5),X1,Y1,'v4'); %real
figure
surf(X1,Y1,Z1,'FaceAlpha',0.3,'FaceColor','r','LineWidth',0.1,'EdgeColor','r')
hold on
surf(X1,Y1,Z2,'FaceAlpha',0.3,'FaceColor','b','LineWidth',0.1,'EdgeColor','b')
%legend('ANN', 'Real');
xlabel('$f$ (Hz)','interpreter','latex') %,'FontSize', 20
ylabel('$a$ (m)','interpreter','latex');
zlabel('$T_c$ ($^{\circ}$C)','interpreter','latex');
view([2.5 10])
grid off
box on
print(gcf,'2.png','-dpng','-r300');

% Time surfaces
x1lin = linspace(min(R1(:,1)),max(R1(:,1)),50);
y1lin = linspace(min(R1(:,2)),max(R1(:,2)),50);
[X1,Y1] = meshgrid(x1lin, y1lin);
Z1 = griddata(R1(:,1),R1(:,2),R2(:,end),X1,Y1,'v4'); %predict
Z2 = griddata(R1(:,1),R1(:,2),R2(:,5),X1,Y1,'v4'); %real
figure
surf(X1,Y1,Z1,'FaceAlpha',0.3,'FaceColor','r','LineWidth',0.1,'EdgeColor','r')
hold on
surf(X1,Y1,Z2,'FaceAlpha',0.3,'FaceColor','b','LineWidth',0.1,'EdgeColor','b')
%legend('RGP', 'Real');
xlabel('$\phi$ (\%)','interpreter','latex');
ylabel('$H_0$ (A/m)','interpreter','latex');
zlabel('$t_\infty$ (s)','interpreter','latex');
view([87.5 10])
grid off
box on
print(gcf,'3.png','-dpng','-r300');

x1lin = linspace(min(R1(:,3)),max(R1(:,3)),50);
y1lin = linspace(min(R1(:,4)),max(R1(:,4)),50);
[X1,Y1] = meshgrid(x1lin, y1lin);
Z1 = griddata(R1(:,3),R1(:,4),R2(:,end),X1,Y1,'v4'); %predict
Z2 = griddata(R1(:,3),R1(:,4),R2(:,5),X1,Y1,'v4'); %real
figure
surf(X1,Y1,Z1,'FaceAlpha',0.3,'FaceColor','r','LineWidth',0.1,'EdgeColor','r')
hold on
surf(X1,Y1,Z2,'FaceAlpha',0.3,'FaceColor','b','LineWidth',0.1,'EdgeColor','b')
%legend('RGP', 'Real');
xlabel('$f$ (Hz)','interpreter','latex');
ylabel('$a$ (m)','interpreter','latex');
zlabel('$t_\infty$ (s)','interpreter','latex');
view([2.5 10])
grid off
box on
print(gcf,'4.png','-dpng','-r300');

