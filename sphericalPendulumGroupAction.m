clear all;
close all;
clc;
        
t0 = 0;
T = 10;
N = 10000;
M = 100;
time = linspace(t0,T,N);
dt = time(2)-time(1);

g = 9.81;
e3 = [0;0;1];
L = 1;

% f = @(input) [input(4:6); g/L * e3]; %output in R6, i.e. the lie algebra of SE(3) with the usual identificatio
f = @(input) [g/L*e3; input(4:6)];
action = @(B,input) [B(1:3,1:3)*input(1:3); B(1:3,1:3)*input(4:6) + hat(B(:,end))*B(1:3,1:3)*input(1:3)];
vecField = @(sigma,p) dexpinvSE3(sigma,f(action(expSE3(sigma),p)));

q0 = rand(3,1);
q0 = q0/norm(q0,2);

v0 = [1;2;0];
v0 = hat(q0)*v0; %so that it leaves in the tangent space at q0
w0 = hat(q0)*v0; 
w = 0.2*w0;

z0 = [q0;w];


field = @(t,z) [hat(z(4:6))*z(1:3); g*hat(e3)*z(1:3)];
[tsol,ysol] = ode45(field, time, z0);

z1 = z0;
z2 = z0;
z3 = z0;

q1 = zeros(3,N);
Energy = zeros(1,N);
q2 = zeros(3,N);
q3 = q2;

AngularMomentum = zeros(1,N);
err = Energy;
err(1) = 0;
omega = q1;

q1(:,1) = q0;
q2(:,1) = q0;
q3(:,1) = q0;
omega(:,1) = w;

En = @(q,w) L^2/2 * norm(hat(w)*q,2)^2 + g*L*e3'*q;
Ang = @(w) L^2*e3'*w;
Energy(1) = L^2/2 * norm(hat(w0)*q0,2)^2 + g*L*e3'*q0;
AngularMomentum(1) = L^2 * e3'*w0;
sigma0 = zeros(6,1);
sigma = sigma0;
it = 1;
for t = time(1:end-1)
    
%     sigma = dt*vecField(sigma0,z); %it is sigma at time dt approximated with Euler method
%     z = action(expSE3(sigma),z);
    
    z1 = LieRK4SPHpendulum(action,f,z1,dt);
    z3 = LieEulerSPHpendulum(vecField,action,z3,dt); %Lie euler
    z2 = RK3SPHpendulum(vecField,action,z1,dt); %3th order Lie RK3
    
    q1(:,it+1) = z1(1:3);
    q2(:,it+1) = z2(1:3);
    q3(:,it+1) = z3(1:3);
    omega(:,it+1) = z1(4:6);
    it = it + 1;    
    err(it) = norm(z1-z2,2);
    Energy(it) = En(q1(:,it),omega(:,it));
    AngularMomentum(it) = Ang(omega(:,it));
    
end


%  %% Solving the system directly on SO(3) and then recovering q(t) on S2
% R = eye(3);
% w = w0;
% qSO3 = zeros(3,N);
% omegaSO3 = qSO3;
% EnergySO3 = zeros(1,N);
% AngularMomentumSO3 = EnergySO3;
% qSO3(:,1) = q0;
% omegaSO3(:,1) = w;
% 
% EnergySO3(1) = L^2/2 * norm(hat(w0)*q0,2)^2 + g*L*e3'*q0;
% AngularMomentumSO3(1) = Ang(w0);
% k = 1;
% 
% %both acting and moved are 3x4 matrices as above
% actionSO3 = @(Acting,Moved) [Acting(1:3,1:3)*Moved(1:3,1:3), Acting(:,end) + Moved(:,end)];
% fSO3 = @(input) [hat(input(:,end)),-g/L * hat(q0) * (input(1:3,1:3)'*e3)]; %matrix of size 3x4, i.e. the Lie algebra of TSO3-> i.e. so(3)xR3
% vecFieldSO3 = @(sigma,p) dexpinvDirect(sigma,fSO3(actionSO3(expDirect(sigma),p)));
% 
% z = [R,w0]; %Matrix of size 3x4, i.e. an element in TSO(3)
% sigma = zeros(6,1);
% 
% for t = time(1:end-1)
%     
%     sigma = dt*vecFieldSO3(sigma,z); %it is sigma at time dt approximated with Euler method
%     z = actionSO3(expDirect(sigma),z);
%     R = z(1:3,1:3);    
%     qSO3(:,k+1) = R * q0; 
%     
%     omegaSO3(:,k+1) = z(:,end);
%     k = k + 1;
%     
%     norm(qSO3(:,k),2)
%     
%     EnergySO3(k) = En(qSO3(:,k),w);
%     AngularMomentumSO3(k) = Ang(omegaSO3(:,k));
%     
% end

% error = vecnorm(q-qSO3);
% figure;
% plot(1:N,error,'r-*')
% title('error against time');

figure;
[x,y,z] = sphere;
% plot3(q1(1,1),q1(2,1),q1(3,1),'r-o',q1(1,end),q1(2,end),q1(3,end),'rd','Markersize',20);
% title('Sol 1')
% hold on;
plot3(q1(1,1:M:end),q1(2,1:M:end),q1(3,1:M:end),'g-*');
% legend('Initial position','Final position','Solution');
% plot3(qSO3(1,1),qSO3(2,1),qSO3(3,1),'g-*','Markersize',20);
% hold on;
% plot3(qSO3(1,:),qSO3(2,:),qSO3(3,:),'r-o');
% hold on;
% plot3(qSO3(1,end),qSO3(2,end),qSO3(3,end),'kd','Markersize',20);

% xlabel('x')
% ylabel('y')
% zlabel('z')
hold on;
% surf(x,y,z)
% colormap('gray');

% figure;
% [x,y,z] = sphere;
% plot3(q2(1,1),q2(2,1),q2(3,1),'k-o',q2(1,end),q2(2,end),q2(3,end),'kd','Markersize',20);
% title('Sol 2')
% hold on;
plot3(q2(1,1:M:end),q2(2,1:M:end),q2(3,1:M:end),'r-d');
hold on;
plot3(q3(1,1:M:end),q3(2,1:M:end),q3(3,1:M:end),'c-o');
hold on;
plot3(ysol(:,1),ysol(:,2),ysol(:,3),'k-','Markersize',5)
legend('RKMK 4', 'RKMK 3', 'Lie Euler','ODE45');
xlabel('x')
ylabel('y')
zlabel('z')

figure;
plot(time,err,'r-*');
title('Difference RKMK 4 and RKMK 3');
figure;
plot(time,Energy,'r-*');
title('Energy RKMK 4');
figure;
plot(time,AngularMomentum,'r-*');
title('Third component of the angular Momentum RKMK 4');

q = ysol(:,1:3);

figure;
nn0 = vecnorm(q');
s0 = string(max(nn0)-min(nn0));
nn1 = vecnorm(q1);
s1 = string(max(nn1)-min(nn1));
nn2 = vecnorm(q2);
s2 = string(max(nn2)-min(nn2));
nn3 = vecnorm(q3);
s3 = string(max(nn3)-min(nn3));

sgtitle("Time norm evolution of the solution")
subplot(2,2,1)
plot(time,nn0,'r-');
title('ODE45, max difference : '+s0);
subplot(2,2,2)
plot(time,nn1,'r-');
title('RKMK 4, max difference : '+s1);
subplot(2,2,3)
plot(time,nn2,'r-');
title('Lie Euler, max difference : '+s2);
subplot(2,2,4)
plot(time,nn3,'r-');
title('RKMK 3, max difference : '+s3);


% figure;
% plot(time,Energy,'r-*',time,EnergySO3,'k-o');
% title('Energy')
% legend('Group action','Standard');
% 
% figure;
% plot(time,AngularMomentum,'r-*',time,AngularMomentumSO3,'k-o');
% title('Angular momentum third component')
% legend('Group action','Standard');

% % Time evolution of the solution
% O = [0;0;0];
% axis(gca,'equal');
% grid on;
% 
% for i = 1:length(time)
%     plot3([O(1), q1(1,i)],[O(2), q1(2,i)],[O(3), q1(3,i)],'r-o');
%     axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]);
%     pause(0.1);
% end