clear all
clc

%in the following line, simple harmonic motion is defined.
theta=0:pi/200:2*pi;
H=0.06; % max rise in [m]
beta=5*pi/6;
gama=beta;
r_r=12; % roller radius in [mm]
r_b=40; %base radius in [mm]
k=20000; %spring constant in [N/m]
delta=10; % spring precompression at s=0 in [mm]

Io2=0.02; %mass moment inertia of link 2 in [kgm^2]
Ibo5=0.05; %mass moment inertia of link 5 in [kgm^2]
a2=0.18; %crank length in [m]
a1=0.3; % distance between two fixed revolute joints
m3=7; % in [kg]
a3=0.1; %in [m]

M=150; % in [Nm]

%defining s
for i=1:length(theta)

if theta(i)<=5*pi/6
    
s(i)=H/2*(1-cos(pi*theta(i)/beta));
else
    s(i)=H-H/2*(1-cos(pi*(theta(i)-gama)/(2*pi-5*pi/6)));
end
end

%defining s'
for i=1:length(theta)
    if theta(i)<=5*pi/6
    s_p(i)=pi*H/(2*beta)*sin(pi*theta(i)/beta);
    else
        s_p(i)=-H/2*(pi/(2*pi-5*pi/6))*sin(pi*(theta(i)-gama)/(2*pi-5*pi/6)); 
    end
end


%defining s''
for i=1:length(theta)
    if theta(i)<=5*pi/6
        s_pp(i)=pi^2*H/(2*beta^2)*cos(pi*theta(i)/beta);
    else
        s_pp(i)=-H/2*(pi/(2*pi-5*pi/6))^2*cos(pi*(theta(i)-gama)/(2*pi-5*pi/6));
    end
end

% figure;
% plot(theta*180/pi, s);
% grid on
% hold on
% plot(theta*180/pi, s_p,"o");
% plot(theta*180/pi, s_pp,"--");
% legend("s","s'", "s''");
% hold off

