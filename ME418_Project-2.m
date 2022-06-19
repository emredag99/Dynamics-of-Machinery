%=====================================================================%
%                    ME418 DoM Cam Mechanism Analysis                 %
%---------------------------------------------------------------------%
%  Coded by:  EMRE DAĞ                                  18/06/2022    %
%=====================================================================%
clear 
close all
clc

%in the following line, simple harmonic motion is defined.
theta=0:pi/200:2*pi;
H=0.06; % max rise in [mm]
beta=5*pi/6; %rise up to beta
gama=beta; %start of return
r_r=12/1000; % roller radius in [m]
r_b=40/1000; %base radius in [m]
k=20000; %spring constant in [N/m]
delta=10/1000; % spring precompression at s=0 in [m]

Io2=0.02; %mass moment inertia of link 2 in [kgm^2]
Ibo5=0.05; %mass moment inertia of link 5 in [kgm^2]
a2=0.18; %crank length in [m]
a1=0.3; % distance between two fixed revolute joints [m]
m3=7; % in [kg]
a3=0.1; %in [m]

g=9.81; % in [m/s^2]
M=150; % during rise period in [Nm] 

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

%plot s, s' and s'' in the same plot
figure;
plot(theta*180/pi, s);
grid on
hold on
plot(theta*180/pi, s_p,"-");
plot(theta*180/pi, s_pp,"--");
legend("s","s'", "s''");
title("Motion of the Cam Follower versus Cam shaft angle [deg]");
hold off

%% check the pressure angle

alpha= atan2(s_p, r_b+r_r+s);
max_alpha=max(alpha);

maxPressureAngle=max(alpha)*180/pi

%plot the pressure angle value
figure
plot(theta*180/pi,alpha*180/pi)
title("Pressure angle versus Cam Shaft Angle")
ylabel("Pressure Angle [deg]")
xlabel("Cam shaft Angle [deg]")
grid on
%% Check undercutting !!!!!

rho= ((r_b+r_r+s).^2+s_p.^2).^(3/2)./((r_b+r_r+s).^2+2*s_p.^2-(r_b+r_r+s).*s_pp)-r_r;
r_follower=ones(length(theta))*r_r;

figure
plot(theta*180/pi, rho)
hold on
plot(theta*180/pi, r_follower)
title("Radius of Curvature and Roller Radius versus Cam Shaft Angle")
ylabel("Radius of Curvature and Follower Radius [m]")
xlabel("Cam Shaft Angle [deg]")
grid on
legend("rho", "r_{follower}")
hold off

%%
phi=acos((a1-a3-r_b-r_r-s)/a2);
g1=s_p./(a2*sin(phi));

g1_p=(a2*s_pp.*sin(phi)-a2*s_p.*cos(phi).*g1)./(a2*sin(phi)).^2;

J=Io2+Ibo5*g1.^2+m3*s_p.^2;
C= Ibo5.*g1.*g1_p+m3*s_p.*s_pp;


for x=1:length(theta)
    if theta(x)<beta
        Q_M(x)=-M.*g1(x);
    else
        Q_M(x)=0;
    end
end
Q_W=-m3*g*s_p;
Q_spr= -k*(s+delta).*s_p;

figure
plot(theta*180/pi, J)
title("Generalized Inertia without Motor and Gearbox versus Generalized Coordinate [deg]")
ylabel("J [kgm^2]")
xlabel("Generalized Coordinate [deg]")
grid on

figure
plot(theta*180/pi, C)
title("Centripetal Inertia without Motor and Gearbox versus Generalized Coordinate [deg]")
ylabel("C [kgm^2]")
xlabel("Generalized Coordinate [deg]")
grid on

figure
plot(theta*180/pi, Q_M)
title("Generalized Force due to External Moment versus Generalized Coordinate [deg]")
ylabel("Q_M [Nm]")
xlabel("Generalized Coordinate [deg]")
grid on

figure
plot(theta*180/pi, Q_W)
title("Generalized Force due to Weight versus Generalized Coordinate [deg]")
ylabel("Q_W [Nm]")
xlabel("Generalized Coordinate [deg]")
grid on

figure
plot(theta*180/pi, Q_spr)
title("Generalized Force due to Spring versus Generalized Coordinate [deg]")
ylabel("Q_spr [Nm]")
xlabel("Generalized Coordinate [deg]")
grid on

slip=1:-0.001:0;
n= 1500*(1-slip); % 


%estimating the maximum torque in the mechanism
%by assuming constant speed at the cam shaft
hiz= 400*pi/30;
T_estimation=C*hiz^2-Q_M-Q_spr-Q_W;
maxTorqueEstimated=max(T_estimation) % estimated maximum torque


%Defining the parameters of the motor
n_n=1421;
T_b=27.27;
s_b=0.315;
a=1.587;
b=1.416;
J_motor=0.0028;
T_motor=T_b./(1+(s_b-slip).^2.*(a./slip-b*slip.^2));

%defining gearbox parameters
r=n_n/(hiz*30/pi); %gear ratio
J_gearbox=(1+0.1*r)*J_motor; %inertia of the gearbox


figure
plot(n, T_motor)
title("Torque-Speed Characteristics of the Motor")
ylabel("T [Nm]")
xlabel("motor speed [rpm]")
grid on

% defining time step and span for Euler method
h=0.00001; %time-step in [s]
t=0:h:3; %span from 0 to 3 seconds

q=zeros(1,length(t));
q_d=zeros(1,length(t));
q_dd=zeros(1,length(t));
J=zeros(1,length(t));
C=zeros(1,length(t));
s=zeros(1,length(t));
s_p=zeros(1,length(t));
s_pp=zeros(1,length(t));
M=zeros(1,length(t));
phi=zeros(1,length(t));
g1=zeros(1,length(t));
g1_p=zeros(1,length(t));

Flywheel=0.06;  % in [kgm^2]
eff=0.9; % gearbox efficiency
counter=1;

for i=1:length(t)-1
    if mod(q(i), 2*pi) <= 5*pi/6
        s(i)=H/2*(1-cos(pi*mod(q(i), 2*pi)/beta));
        s_p(i)=pi*H/(2*beta)*sin(pi*mod(q(i), 2*pi)/beta);
        s_pp(i)=pi^2*H/(2*beta^2)*cos(pi*mod(q(i),2*pi)/beta);
        M(i)=150;
       
    else
        
        s(i)=H-H/2*(1-cos(pi*(mod(q(i), 2*pi)-gama)/(2*pi-5*pi/6)));
        s_p(i)=-H/2*(pi/(2*pi-5*pi/6))*sin(pi*(mod(q(i), 2*pi)-gama)/(2*pi-5*pi/6));
        s_pp(i)=-H/2*(pi/(2*pi-5*pi/6))^2*cos(pi*(mod(q(i), 2*pi)-gama)/(2*pi-5*pi/6));
        M(i)=0;
    end
    
    
    phi(i)=acos((a1-a3-r_b-r_r-s(i))/a2);
    g1(i)=s_p(i)./(a2*sin(phi(i)));
    g1_p(i)=(a2*s_pp(i).*sin(phi(i))-a2*s_p(i).*cos(phi(i)).*g1(i))./(a2*sin(phi(i))).^2;
    J(i)=Io2+Ibo5*g1(i).^2+m3*s_p(i).^2+J_gearbox+J_motor*r^2+Flywheel*r^2;
    C(i)=Ibo5.*g1(i).*g1_p(i)+m3*s_p(i).*s_pp(i);
    Q_M=-M(i).*g1(i);
    Q_W=-m3*g*s_p(i);
    Q_spr=-k*(s(i)+delta).*s_p(i);
    slip=1-abs(q_d(i))*r/(1500*pi/30);
    Q_motor(i)= T_b./(1+(s_b-slip).^2.*(a./slip-b*slip.^2))*r*eff;
    Q(i)=Q_motor(i)+Q_M+Q_W+Q_spr;
    
    alpha(i)=atan2(s_p(i),r_b+r_r+s(i));
    
    q(i+1)=q(i)+h*q_d(i);
    q_dd(i)= (Q(i)-C(i)*q_d(i)^2)/J(i);
    q_d(i+1)= q_d(i)+h*q_dd(i);
    
    if q_d(i)>400*pi/30 & counter==1
        steadtStateTime=t(i);
        counter=71;   
    end
end

%% creating expression for fluctuation
q_dd_avg=400*pi/30;
fluctuation=(max(q_d(ceil(length(t)/3):length(t)-1))-min(q_d(ceil(length(t)/3):length(t)-1)))/q_dd_avg*100;
%% creating expression for stability
omega_bd=(1-s_b)*max(n); 
omega_min=min(q_d(ceil(length(t)/3):length(t)-1))*30/pi; %[rpm]
omega_max=max(q_d(ceil(length(t)/3):length(t)-1))*30/pi; %[rpm]
stability= (r*omega_min-omega_bd)/(max(n)-omega_bd)*100;
%% running cost
omega_Avg=(omega_min+omega_max)/2;
runningCost= abs(omega_Avg*r-n_n)/(max(n)-omega_bd)*100;
%% performance parameters
steadtStateTime,stability, fluctuation,runningCost

%%
figure
plot(t,abs(q_d)*30/pi)
title("Cam Shaft Speed versus Time")
ylabel("Cam Shaft Speed [rpm]")
xlabel("Time [s]")
grid on

figure
plot(q_d(1:length(Q_motor))*30/pi, Q_motor)
title("Motor Torque Reduced to the Cam Shaft versus Its Speed")
ylabel("Motor Torque [Nm]")
xlabel("Cam Shaft Speed [rpm]")
grid on

%plotting the torque due to generalized inertia versus time
figure
plot(t, J.*q_dd,"r")
title("Torque due Generalized Inertia versus Time")
ylabel("Torque due Generalized Inertia [Nm]")
xlabel("Time [s]")
grid on

%plotting torque due to centripetal inertia
figure
plot(t, C.*q_d.^2,"k")
title("Torque due to Centripetal Inertia versus Time")
ylabel("Torque due to Centripetal Inertia [Nm]")
xlabel("Time [s]")
grid on

%% contact force
alpha(length(alpha)+1)=alpha(length(alpha));
F45y=(Ibo5*(g1.*q_dd+g1_p.*q_d.^2)+M)./(a2*sin(phi));
F_contact=(m3*(s_p.*q_dd+s_pp.*q_d.^2)+m3*g+k*(s+delta)+F45y)./(cos(alpha));

minimumContactForce=min(F_contact)

figure
plot(t, F_contact,"k")
title("Contact Force versus Time")
ylabel("Contact Force [N]")
xlabel("Time [s]")
grid on

