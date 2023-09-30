%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLANE WAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

%info
version="Synch";
fnt="Latin Modern Mono Slanted";
mrksz=8;
fntsz=14;

% laser
a0=1;
w0=1;

%particle
gm0=5; % initial gamma
beta=sqrt(1-1./gm0.^2);
v0=[beta(1),0,0]';

wi=0.001;
wf=5*w0+0.001;
Nw=250;
w=linspace(wi,wf,Nw); % frequency array

Thdim=1;
Thmax=0; % in radians

[mat,t,r]=zfunc(v0,a0,w0,wi,wf,Nw,Thmax,Thdim);
toc

%PLOTS
tic
hold on;

%theta=0;
output_precision(2);

%1D
theta=0; %suppose
%plot(w/w0,mat(1:end,floor(Thdim/2))/max(mat(1:end,floor(Thdim/2))),'b',"markersize", mrksz);%(1:end,floor(Thdim/2))/max(mat(1:end,floor(Thdim/2)));
plot(w/w0,mat(1:end,1)/max(mat(1:end,1)),'b',"markersize", mrksz);%(1:end,floor(Thdim/2))/max(mat(1:end,floor(Thdim/2)));
xlabel("$\\omega/\\omega_0$",'rotation',0,"interpreter", "latex",'fontsize',fntsz);
ylabel("$\\frac{d^2I}{d\\omega d\\Omega} (\\theta=0)$",'rotation',90,"interpreter", "latex",'fontsize',fntsz);
title(["$" sprintf("a_0=%u, \\omega_0=%u,\\gamma_0=%u",a0,w0,gm0) "$"],"interpreter", "latex",'fontsize',fntsz*1.7);
hold on;
xlim([0 5])

figure();
plot(t,r(1,:),'b',"markersize", mrksz);
xlabel("$t[1/\\omega_0]$",'rotation',0,"interpreter", "latex",'fontsize',fntsz);
ylabel("$x[c/\\omega_0]$",'rotation',90,"interpreter", "latex",'fontsize',fntsz);
title(["$" sprintf("a_0=%u, \\omega_0=%u,\\gamma_0=%u",a0,w0,gm0) "$"],"interpreter", "latex",'fontsize',fntsz*1.7);

figure();
plot(r(1,:),r(2,:),'b',"markersize", mrksz);
xlabel("$x[c/\\omega_0]$",'rotation',0,"interpreter", "latex",'fontsize',fntsz);
ylabel("$y[c/\\omega_0]$",'rotation',90,"interpreter", "latex",'fontsize',fntsz);
title(["$" sprintf("a_0=%u, \\omega_0=%u,\\gamma_0=%u",a0,w0,gm0) "$"],"interpreter", "latex",'fontsize',fntsz*1.7);

#{
figure();
plot3(r(1,:),r(2,:),r(3,:),'b',"markersize", mrksz);
xlabel("$x[c/\\omega_0]$",'rotation',0,"interpreter", "latex",'fontsize',fntsz);
ylabel("$y[c/\\omega_0]$",'rotation',90,"interpreter", "latex",'fontsize',fntsz);
zlabel("$z[c/\\omega_0]$",'rotation',90,"interpreter", "latex",'fontsize',fntsz);
#}

%2D
figure();
colormap ("default");
colormap ("jet");
plt1=pcolor (mat); 
shading interp;
axis ("tic", "square","labelx","tight","labely","tight");

xlabel("$\\gamma \\theta_{detec} (rad)$",'rotation',0,"interpreter", "latex");
set(gca, 'xTick', [1, floor(Thdim/2), Thdim-1]);
set(gca, 'xTickLabel', [-pi, 0, pi]);

ylabel("$\\omega/\\omega_c$",'rotation',90,"interpreter", "latex");
set(gca, 'yTick', linspace(1,Nw,9));
set(gca, 'yTickLabel', linspace(0,wf/w0,9));

title(["$\\frac{d^2I}{d\\omega d\\Omega}(" sprintf("a_0=%u, \\gamma_0=%u)",a0,gm0) "$"] ,'rotation',0,"interpreter", "latex",'fontsize',fntsz*1.7);

#}
  #{
  hold on;
  positionVector2 = [0.1, 0.1, 0.4, 0.2];    % position of second subplot
  subplot('Position',positionVector2);
  plot (linspace(mmin,mmax,mdim), vecO/max(vecO), 'b',"markersize", mrksz);
  xlabel("$\\omega/\\omega_0$",'rotation',0,"interpreter", "latex");
  ylabel("$dI/d\\omega$",'rotation',90,"interpreter", "latex");
  axis("tight");

  hold on;
  positionVector3 = [0.501, 0.3,0.2,0.5];    % position of second subplot
  subplot('Position',positionVector3);
  plt3=plot (linspace(Thmin,Thmax,Thdim), fliplr(vecdI'/max(vecdI)), 'b',"markersize", mrksz);
  direction = [0 0 1];
  rotate(plt3,direction,-90);
  axis ("tic", "labelx","tight");
  lbl=xlabel("$dI/d\\Omega$",'rotation',0,"interpreter", "latex");
  #}
  
  #{
#printing
print -depslatexstandalone sombrero;
## process generated files with pdflatex
system ("latex sombrero.tex");
## dvi to ps
system ("dvips sombrero.dvi");
## convert to png
system ("gs -dNOPAUSE -dBATCH -dSAFER -sDEVICE=png16m -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -r100x100 -dEPSCrop -sOutputFile=sombrero.png sombrero.ps");
system (["cp sombrero.png sync.png"]);
system ("rm sombrero.aux sombrero.dvi sombrero.log sombrero.ps sombrero.tex sombrero-inc.eps sombrero.png");
system ("clear");
#}
