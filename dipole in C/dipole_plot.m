close
clear

filename=('dipole.csv');
M = csvread(filename);

figure
plot3(M(:,2),M(:,3),M(:,4),'Linewidth',2); hold on; 
set(gca,'TickLabelInterpreter','latex','Fontsize',14)

ylabel('$ y / d_{\rm p} $','Interpreter','latex','Fontsize',20);
xlabel('$ x / d_{\rm p} $','Interpreter','latex','Fontsize',20);
zlabel('$ z / d_{\rm p} $','Interpreter','latex','Fontsize',20);

