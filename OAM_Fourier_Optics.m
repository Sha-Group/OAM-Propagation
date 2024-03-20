
clc;clear

%% configuration
wavelength=1;  %  ***wavelength
k0=2*pi/wavelength;  %  wavenumber
Z=120*pi;  % wave impedance
Volume=1;
%an=(Volume/(pi*4/3)).^(1/3);
TL=2;   % Topological charge   ***

%% number and position of source points
N=12;   %  *** number of point sources

% radius and angle
radius= 0.5;  %  *** radius
theta=linspace(0,2*pi-2*pi/(N),N);
[x00,y00]=pol2cart(theta,radius);

%% sampling the observation area xoy plane
%  z direction slices
delta=0.4;       % step ***
size=30;         % xoy range ***

% xoy plane
xo=-size:delta:size-delta;
yo=-size:delta:size-delta;
zo=1;  % observation z  ***

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
xlabel('x (a.u.)')
ylabel('y (a.u.)')
axis equal
% N_p=index;
title('Location of point sources')

%% dyadic Green function
Greenxx=zeros(N_p,N);
Greenyy=zeros(N_p,N);
Greenzz=zeros(N_p,N);
Greenxy=zeros(N_p,N);
Greenxz=zeros(N_p,N);
Greenyz=zeros(N_p,N);

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
        const1=-j*k0*Z*k0*Volume*exp(-j*alpha)/(4*pi*alpha^3);
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

%  dyadic Green function
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

%  transverse component
[x,y]=meshgrid(xo,yo);
[phi, rho]=cart2pol(x,y);
Ex_re=reshape(Ex,Ny,Nx);
%Ex_re=Ex_re./norm(Ex_re,'fro');

figure(2)

subplot(2,1,1)
pcolor(x,y,real(Ex_re));
shading interp;
colorbar;
axis equal
axis([- size size-delta -size size-delta])
title('Re(E_x) at an initial plane')

subplot(2,1,2)
pcolor(x,y,abs(Ex_re));
shading interp;
colorbar;
axis equal
axis([- size size-delta -size size-delta])
title('|E_x| at an initial plane')

%% 

% initial parameters for Fourier optics
%lambda = 632.8e-9; % wavelength
%k = 2*pi / lambda; % wavenumber
Scale=5;            % scaling factor for zero padding
L = Scale*2*size;   % total length
M = Scale*Nx;       % total points
z = 5;              % propagating distance *** 

% 网格定义
dx = L / M;
% x = -L/2 : dx : L/2 - dx;
% [X, Y] = meshgrid(x, x);

% initial vortex beam (in real space)
U0 = Ex_re;
Upad=padarray(U0,[(M-Nx)/2,(M-Nx)/2],0,'both'); % zero padding

% angular spectrum representation of vortex beam
F_U0 = fftshift(fft2((Upad)));

% angular frequency coordinates
df = 1 / (L);
fX = -1/(2*dx) : df : 1/(2*dx) - df;
[ffX, ffY] = meshgrid(fX, fX);

% angular propagation
flag=(1 - wavelength^2 * (ffX.^2 + ffY.^2))>=0;
%sum(sum(flag))
H = exp(-j * k0* z * conj(sqrt(1 - wavelength^2 * (ffX.^2 + ffY.^2))));
F_Uz = H .* F_U0;

% inverse Fourier transform
Uz = (ifft2(fftshift(F_Uz)));
Uz = Uz((M-Nx)/2+1:(M+Nx)/2,(M-Nx)/2+1:(M+Nx)/2);

%% Solution by Green's function method
%  z direction slices
%delta=0.2;      % step ***
%size=8;         % xoy range ***

% xoy plane
xo=-size:delta:size-delta;
yo=-size:delta:size-delta;
zo=zo+z;  % observation z

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
        const1=-j*k0*Z*k0*Volume*exp(-j*alpha)/(4*pi*alpha^3);
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
    Jx(index_s)=exp(i*TL*theta(index_s));
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
%Ex_re=Ex_re./norm(Ex_re,'fro');

figure(3)

subplot(2,1,1)
pcolor(x,y,abs(Ex_re));
shading interp;
colorbar;
axis equal
axis([- size size-delta -size size-delta])
title('|E_x| by Green fuction method')

subplot(2,1,2)
pcolor(x,y,abs(Uz));
%rpcolor(abs(Uz));
shading interp;
colorbar;
axis equal
axis([- size size-delta -size size-delta])
title('|E_x| by Fourier optics')

% calculation error
norm(Uz-Ex_re,'fro')/norm(Ex_re,'fro')
