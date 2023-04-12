%Finite Elements Term Project
%Anthony Thompson
%Jude Anike

clear
clc
E = 2.5*10^10;          %Modulus of elasticity
rho=2500;               %beam density kg/m^3
A = 0.25*0.4;           %Cross section area
I = 0.25*(0.4^3)/12;    %Area moment of inertia
G = E/2.5;              %Shear Modulus
ks = 5/6;               %Timeshenko shear coefficient
dt = 0.001;             %size of time step
Ntimestep = 1000;       %number of time steps
L = 10;                 %length of the beam
n = 100;                %number of elements
h=L/n;                  %Size of each element
f=2;                    %Applied load frequency in Hz
t=0;
load=100;               %Amplitude of the distributed load
loadtime=1000;            %duration load is applied in milliseconds
size=(load/100);
%%
syms ksi;               
phi1=-(9/16)*(ksi+(1/3))*(ksi-(1/3))*(ksi-1);   %Cubic shape function 1
phi2=(27/16)*(ksi+1)*(ksi-(1/3))*(ksi-1);       %Cubic shape function 2
phi3=-(27/16)*(ksi+1)*(ksi+(1/3))*(ksi-1);      %Cubic shape function 3
phi4=(9/16)*(ksi+1)*(ksi+(1/3))*(ksi-(1/3));    %Cubic shape function 4
phi=[phi1 phi2 phi3 phi4];                          %vector for calling cubic shape functions
dphi=diff(phi,ksi);                       %vector for calling derivative of cubic shape functions

M1E=[];                                         %element Mass Matrix 1
for i=1:4                               %Loop that creates element Mass Matrix 1
    for j=1:4
        M1E(i,j)=(h/2)*rho*A*int(phi(i)*phi(j),-1,1);
    end
end

M1G=zeros(3*n+1);                       %Creating the cubic part of the Global Mass Matrix
for i=1:3:3*n-2
    M1G(i:i+3,i:i+3)=M1G(i:i+3,i:i+3)+M1E;
end
             
g1=(1/2)*ksi*(ksi-1);                     %Quadratic shape function 1
g2=(1-ksi^2);                              %Quadratic shape function 2
g3=(1/2)*ksi*(ksi+1);                     %Quadratic shape function 3
g=[g1 g2 g3];                               %vector for calling quadratic shape functions
dg=diff(g);                         %Vector for calling derivative of quadratic shape functions

M2E=[];                                        %Element Mass Matrix 2
for i=1:3                               %Loop that creates element Mass Matrix 2
    for j=1:3
        M2E(i,j)=(h/2)*rho*I*int(g(i)*g(j),-1,1);
    end
end

M2G=zeros(2*n+1);                       %Creating the quadratic part of the Global Mass Matrix
for i=1:2:2*n-1
    M2G(i:i+2,i:i+2)=M2G(i:i+2,i:i+2)+M2E;
end

Mst=zeros(5*n+2);                       %Creating the Mass state Matrix
Mst(1:3*n+1,1:3*n+1)=M1G;
Mst(3*n+2:5*n+2,3*n+2:5*n+2)=M2G;
Mst=sparse(Mst);

K1E=[];                        %Creating top left quadrant of element stiffness matrix
for i=1:4
    for j=1:4
     K1E(i,j)=(2/h)*G*A*ks*int(dphi(i)*dphi(j),-1,1);  
    end
end

K2E=[];                        %Creating top right quadrant of element stiffness matrix
for i=1:4
    for j=1:3
     K2E(i,j)=-G*A*ks*int(dphi(i)*g(j),-1,1);  
    end
end 

K3E=[];                        %Creating bottom left quadrant of element stiffness matrix
for i=1:3
    for j=1:4
     K3E(i,j)=-G*A*ks*int(g(i)*dphi(j),-1,1);  
    end
end

K4E=[];                         %Creating bottom right quadrant of element stiffness matrix
K4Ea=[];                        
K4Eb=[];

for i=1:3
    for j=1:3
     K4Ea(i,j)=E*I*(2/h)*int(dg(i)*dg(j),-1,1);
     K4Eb(i,j)=G*A*ks*(h/2)*int(g(i)*g(j),-1,1);
     K4E(i,j)=K4Ea(i,j)+K4Eb(i,j);
    end
end 

K1G=zeros(3*n+1);                       %Creating the top left quadrant of Global stiffness Matrix
for i=1:3:3*n-2
    K1G(i:i+3,i:i+3)=K1G(i:i+3,i:i+3)+K1E;
end

K2G=zeros(3*n+1,2*n+1);                       %Creating the top right quadrant of Global stiffness Matrix
j=1;
for i=1:3:3*n-2
    K2G(i:i+3,j:j+2)=K2G(i:i+3,j:j+2)+K2E;
    j=j+2;
end

K3G=zeros(2*n+1,3*n+1);                       %Creating the bottom left quadrant of Global stiffness Matrix
j=1;
for i=1:2:2*n-1
    K3G(i:i+2,j:j+3)=K3G(i:i+2,j:j+3)+K3E;
    j=j+3;
end

K4G=zeros(2*n+1,2*n+1);                       %Creating the bottom right quadrant of Global stiffness Matrix
for i=1:2:2*n-1
    K4G(i:i+2,i:i+2)=K4G(i:i+2,i:i+2)+K4E;
end

Kst=zeros(5*n+2,5*n+2);                         %Combining global stiffness matrix
Kst(1:3*n+1,1:3*n+1)=K1G;
Kst(1:3*n+1,3*n+2:5*n+2)=K2G;
Kst(3*n+2:5*n+2,1:3*n+1)=K3G;
Kst(3*n+2:5*n+2,3*n+2:5*n+2)=K4G;

Kst(:,3*n+1)=[];                        %Applying boundary conditions
Kst(3*n+1,:)=[];
Kst(:,1)=[];
Kst(1,:)=[];
Mst(:,3*n+1)=[];
Mst(3*n+1,:)=[];
Mst(:,1)=[];
Mst(1,:)=[];

Mst=sparse(Mst);
Kst=sparse(Kst);

q=-load*sin(2*pi*f*t);                       %Creating Force Element Vector
ForceE=[];
for i=1:4                                 
    ForceE(i)=int(phi(i),-1,1)*(h/2);
end
ForceE=ForceE';                             

Qst=zeros(5*n+2,1);                           %Creating the full Force State vector
for i=1:3:3*n-2
    Qst(i:i+3)=Qst(i:i+3)+ForceE;
end
Qst(3*n+1,:)=[];                    %Applying boundary conditions
Qst(1,:)=[];
Qst=sparse(Qst);

ui=zeros(5*n,1);                    %setting initial conditions ui=initial position state vector
dui=zeros(5*n,1);                   %dui=initial velocity state vector
ddui=inv(Mst)*(q*Qst);              %ddui=initial acceleration state vector

u=[];
du=[];
ddu=[];

U=[];                           %State vector that will be used for plotting


for i=1:Ntimestep
    t=i*dt;
    if i<loadtime                   %setting duration of the loading
    q=-load*sin(2*pi*f*t);
    else
        q=0;
    end
    COMP=[Mst+(Kst*(dt^2)*(1/4))];
    ddu=inv(COMP)*((q*Qst)-(Kst*(ui+dui*dt+(1/4)*ddui*dt^2)));
    du=dui+(.5*ddui+.5*ddu)*dt;
    u=ui+dui*dt+.5*(.5*ddui+.5*ddu)*dt^2;
    
    
    
    
    ui=u;           %setting state solution vector equal to initial for next time step
    dui=du;
    ddui=ddu;
     mid_point_soln_arr(i) = u((3*n)/2-1,1);
     U(:,:,i)=u;                %3D matrix that holds displacement vectors for each time step
end

time_arr = [dt:dt:dt*Ntimestep]; 
plot(time_arr, mid_point_soln_arr,'ro');

fixednode=zeros(1,1,Ntimestep);     
U1=U(1:299,:,:);                %cutting just displacement vector from full state vector
U1=[fixednode;U1;fixednode];    %putting fixed nodes back in displacement vector

x=0:3*n;        
xh=x*(1/3)*h;                   %vector for plot



figure
title([num2str(load), ' Newton distributed load'])
xlabel('horizontal position in meters')
ylabel('displacement in meters')
axis([0 10 -6.2*10^-4*size 6*10^-4*size])
zz=0;
while zz==0;
demo=input('type y and enter when ready to run demo. ','s');
switch(demo)
    case 'y'
for i=1:Ntimestep
 
   clf
   title([num2str(load),' Newton distributed load.  t=',num2str(i/1000),' seconds'])
xlabel('horizontal position in meters')
ylabel('displacement in meters')
   axis([0 10 -6.2*10^-4*size 6*10^-4*size])
    axis manual
    hold on
    plot(xh,U1(:,:,i))
    %daspect([1 1 1]);
    pause(1/120)
    hold off
   
 
end
zz=1;
again=input('would you like to run the demo again?  [y,n] ','s');
switch(again)
    case 'y'
        zz=0;
end
        
end
end

% for i=1:Ntimestep
%  
%    clf
%    axis([0 10 -10*10^-4 10*10^-4])
%     axis manual
%     hold on
%     plot(xh,U1(:,:,i))
%     %daspect([1 1 1]);
%     pause(1/120)
%     hold off
%    
%  
% end




                


