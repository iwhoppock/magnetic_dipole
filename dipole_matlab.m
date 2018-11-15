
%This models a charged particle in a magnetic mirror/magnetic dipole. For
%it to work properly, the initial magnetid field at the starting location
%of the particle must be 1. This may be seen if the user uncomments the if
%statement in the time loop. Changes to the initial position will likely
%change the parameters of the magnetic field (all one must do is print out
%the initial field, and then multiply/divide the components of said field
%to make it one.) The units are normalised by the: speed is in terms of the
%Alfven spee, frequencies the proton cyclotron frequency, lengths the
%proton inertial length, etc. To change to an alpha, replace mass with 4
%and charge with 2. Change the velocities gradually. This code implements
%the Boris algorithm for non-relativistic particle tracing. The quantity
%vAc is the ratio of the alfven speed by the speed of light (c). 

%NB: this is not a fast algorithm because Matlab. For longer runs, use the
%code written in C. 

clear
close

dt = 0.01;
mass = 1.0;
charge = 1.0;
vAc = 3e-4;

duration = 200000;

X = zeros(duration,3);
V = zeros(duration,3);


v = [0, 0.1, 0.1];
x = [10, 0, 0];

E = [0., 0., 0.];

for time = 1:1:duration
    B = dipole(x);
%     if time == 1
%         B
%     end
    [x,v] = boris_rotaion(x,v,charge,mass,vAc,dt,B,E);
    X(time,:) = x;
    V(time,:) = v;
end    

plot3(X(:,1),X(:,2),X(:,3),'k','Linewidth',2); hold on;
set(gca,'TickLabelInterpreter','latex','Fontsize',14)

ylabel('$ y / d_{\rm p} $','Interpreter','latex','Fontsize',20);
xlabel('$ x / d_{\rm p} $','Interpreter','latex','Fontsize',20);
zlabel('$ z / d_{\rm p} $','Interpreter','latex','Fontsize',20);

function [x,v] = boris_rotaion(x,v,charge,mass,vAc,dt,B,E)
    t = charge ./ mass .* B .* 0.5 .* dt;
    s = 2 .* t ./ (1 + t.*t);
    v_minus = v + charge ./ (mass * vAc) .* E .* 0.5 .* dt; 
    v_prime = v_minus + cross(v_minus,t);
    v_plus = v_minus + cross(v_prime,s);
    v = v_plus + charge ./ (mass * vAc) .* E .* 0.5 .* dt;
    x = x + v .* dt;
end

function B = dipole(x)
    position_mag = sqrt(dot(x,x));
    m = 1000.;
    B(1) = (3 * m * x(1) * x(3)) / (position_mag^5);
    B(2) = (3 * m * x(2) * x(3)) / (position_mag^5);
    B(3) = ((m) / (position_mag^3)) * (( (3 * x(3) * x(3)) / (position_mag^2) ) - 1);
end