clc;
close all;
clear  all;

%% configuration
wavelength=1;  %  ***wavelength
k0=2*pi/wavelength; % wavenumber
Z=120*pi; % wave impedance
Volume=1;
an=(Volume/(pi*4/3)).^(1/3);
TL=3;   % Topological charge

%% number and position of source points
N=12;   %  *** number

% radius and angle
radius= 0.5;  %  *** radius
theta=linspace(0,2*pi-2*pi/(N),N);
[x00,y00]=pol2cart(theta,radius);

%% sampling the observation area at xoy plane

delta=0.2;      % step ***
size=8;         % xoy range ***

% xoy plane
xo=-size:delta:size-delta;
yo=-size:delta:size-delta;
zo=2;  % z coordinate of observation plane

% size
Nx=length(xo);
Ny=length(yo);

%  position of observation point
index=0;
for m=1:Nx;  %  X
    for n=1:Ny;  %  Y
            index=index+1;
            posx(index)=xo(m);
            posy(index)=yo(n);
            posz(index)=zo(1);
    end
end
N_p=index;

figure(1);
plot(x00,y00,'r.')
hold on
plot([x00,0.5],[y00,0],'b--')
axis equal
% N_p=index;
% title('observation point')
title('Location of point sources')

%% dyadic Green function
for index_f=1:N_p; % field point number
    %  counting
    %N_p-index_f
    
    %  position of observation point
    x=posx(index_f);
    y=posy(index_f);
    z=posz(index_f);
    
    for index_s=1:N;
        
        %  position of source point
        xx=x00(index_s);
        yy=y00(index_s);
        zz=0;
        
        %  distance
        R=sqrt((xx-x)^2+(yy-y)^2+(zz-z)^2);
        alpha=k0*R;
        %  direction
        cosx=(x-xx)/R;
        cosy=(y-yy)/R;
        cosz=(z-zz)/R;
        %  costant
        const1=j*k0*Z*k0*Volume*exp(-j*alpha)/(4*pi*alpha^3);
        const2=3-alpha^2+3*j*alpha;
        const3=(alpha)^2-1-j*alpha;
        Greenxx(index_f,index_s)=const1*(const3+cosx*cosx*const2);
        Greenyy(index_f,index_s)=const1*(const3+cosy*cosy*const2);
        Greenzz(index_f,index_s)=const1*(const3+cosz*cosz*const2);
        Greenxy(index_f,index_s)=const1*cosx*cosy*const2;
        Greenxz(index_f,index_s)=const1*cosx*cosz*const2;
        Greenyz(index_f,index_s)=const1*cosy*cosz*const2;
        
    end
end

%  dyadic green function
G=[Greenxx Greenxy Greenxz;
    Greenxy Greenyy Greenyz;
    Greenxz Greenyz Greenzz];

%  decomposition of source point (x polaried OAM source)
for index_s=1:N;
    Jx(index_s)=exp(i*TL*theta(index_s));  % helical wavefront
    Jy(index_s)=0;
    Jz(index_s)=0;
end

J=[Jx Jy Jz].';

% radiated field
Etot=G*J;

Ex=Etot([1:N_p],1);
Ey=Etot([N_p+1:2*N_p],1);
Ez=Etot([2*N_p+1:3*N_p],1);

%  transverse component  ***
[x,y]=meshgrid(xo,yo);
[phi, rho]=cart2pol(x,y);
Ex_re=reshape(Ex,Ny,Nx);
Ex_re=Ex_re./norm(Ex_re,'fro');

%  plot field on transverse plane
figure(2);

subplot(2,2,1)
pcolor(x,y,abs(Ex_re));
shading interp;
colorbar;
axis equal
axis([- size size-delta -size size-delta])
title('|E_x|')

subplot(2,2,2)
pcolor(x,y,real(Ex_re));
shading interp;
colorbar;
axis equal
axis([- size size-delta -size size-delta])
title('Re(E_x)')

%% derivative with x
z=fft(Ex_re);      % column fft
z=fftshift(z,1);   % column shift

n=0:Nx-1;
v=j*pi/delta*(n-(Nx)/2)/((Nx)/2); % j*x in fourier space
z=z.*v.';  % column
z=ifftshift(z,1);  % column inverse shift

y_x=(ifft(z)); % column ifft


%% derivative with y
z=fft(Ex_re.');      % transpose and column fft
z=fftshift(z,1);     % column shift

n=0:Ny-1;
v=j*pi/delta*(n-(Ny)/2)/((Ny)/2); % j*y in fourier space
z=z.*v.';  % column
z=ifftshift(z,1);  % column inverse shift
y_y=(ifft(z)); % column ifft

y_y=y_y.'; % transpose return back

%% ∂‘ Phi«Ûµº

y_phi=y_x.*x-y_y.*y;
p_den=conj(Ex_re).*(y_phi)*(-i);

%figure(3)
subplot(2,2,3)
pcolor(x,y,real(p_den));
shading interp;
colorbar;
axis equal
axis([- size size-delta -size size-delta])
title('Re(L_z^f)')

subplot(2,2,4)
pcolor(x,y,imag(p_den));
shading interp;
colorbar;
axis equal
axis([- size size-delta -size size-delta])
title('Im(L_z^f)')

sum(sum(p_den))
