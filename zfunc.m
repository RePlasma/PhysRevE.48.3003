%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLANE WAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FWM,t,r]=zfunc(v0,a0,w0,wi,wf,Nw,Thmax,Thdim)
% initial definitions
q=1; m=1; c=1;                                                                        
dt=0.01; % 0.1 is probably too coarse
No=12; % #total simulation periods
Nenv=No-1; % #envelope periods
gm0=1/sqrt(1-v0'*v0/c^2);
vabs=sqrt(v0'*v0);
ti=0; tf=round(2*pi*No); % because q=m=B=c=1

% initial arrays
t=(ti+dt):dt:tf;
w=linspace(wi,wf,Nw);
N=round((tf-ti)/dt); % #time steps

% allocation
u=zeros(3,N);
v=zeros(3,N);
v(:,1)=v0;
r=zeros(3,N);

% initial conditions
u(:,1)=(v0./sqrt(1-v0'*v0/c^2));
r(:,1)=[0,0,0]';

% wave travelling in the -x direction
% Envelope
nenv=Nenv*round(2*pi/dt);
for n=2:nenv
    f = sin(n*dt/(Nenv*4))^2; % envelope function
    E0 = f*[0,a0*w0*sin(n*dt+r(1,n-1)),0]'; % Electric field
    B = f*[0,0,a0*w0*sin(n*dt+r(1,n-1))]'; % Magnetic field
    um = u(1:3,n-1)+E0*q*dt/(2*m);
    gm = sqrt(1+um'*um);
    co = q*dt/(m*gm);
    up = um+cross(um,B)*co/2;
    uP = um+(co/(1+(co/2)^2))*cross(up,B);
    u(1:3,n)=uP+q*dt*E0/(2*m);
    v(:,n)=(c*u(:,n)./sqrt(c^2+u(:,n)'*u(:,n)));
    r(1:3,n)=r(1:3,n-1)+dt*v(1:3,n);
end
% Boris pusher algorithm
for n=nenv+1:N
    E0 = [0,a0*w0*sin(n*dt+r(1,n-1)),0]'; % Electric field
    B = [0,0,a0*w0*sin(n*dt+r(1,n-1))]'; % Magnetic field
    um = u(1:3,n-1)+E0*q*dt/(2*m);
    gm = sqrt(1+um'*um);
    co = q*dt/(m*gm);
    up = um+cross(um,B)*co/2;
    uP = um+(co/(1+(co/2)^2))*cross(up,B);
    u(1:3,n)=uP+q*dt*E0/(2*m);
    v(:,n)=(c*u(:,n)./sqrt(c^2+u(:,n)'*u(:,n)));
    r(1:3,n)=r(1:3,n-1)+dt*v(1:3,n);
end
beta = v/c;   %Beta
pbeta = zeros(3,N);
for k=1:N-1
    pbeta(:,k) = (beta(:,k+1)-beta(:,k))/dt;  % beta derivative
end
pbeta(:,N) = pbeta(:,N-1);   % to keep array size N 

LeftIntegrand = zeros(3,N); %
RightExp = zeros(1,N);
FW = zeros(1,Nw);
FWM = zeros(Nw,Thdim);

for nth=1:Thdim
    th = -Thmax+2*(nth-1)*Thmax/Thdim;%
    parfor k=nenv:N
        nv = [cos(th); sin(th); 0];
        LeftIntegrand(:,k) = cross(nv,cross(nv-beta(:,k),pbeta(:,k)))./(1-nv'*beta(:,k))^2;
        RightExp(k) = 1i*(t(k)-nv'*r(:,k)/c);
    end
    parfor l=1:Nw
        AllIntegrand = bsxfun(@times,LeftIntegrand,exp(w(l)*RightExp));
        Integral = sum(AllIntegrand,2)*dt;
        FW(l) = real((Integral'*Integral)); %
    end
    FWM(1:end,nth)=FW';
endfor

