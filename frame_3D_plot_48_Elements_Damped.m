 %Anthony Thompson

clear
clc
E = 2.5*10^10;          %Modulus of elasticity
rho=2500;               %beam density kg/m^3
A = 0.25*0.4;           %Cross section area
I = 0.25*(0.4^3)/12;    %Area moment of inertia
G = E/2.5;              %Shear Modulus
ks = 5/6;               %Timeshenko shear coefficient
mu=4000;                %Viscous Damping factor
dt = 0.001;             %size of time step
Ntimestep = 1500;       %number of time steps
L = 10;                 %length of linearly connected elements
n = 40;                 %number of linearly connected elements
ncm=8;                  %number of cross member elements
h=L/n;                  %Size of each element
f=2;                    %Applied load frequency in Hz
t=0;
load=250;               %Amplitude of the distributed load
loadtime=150;            %duration load is applied in milliseconds
size=1;
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



M1G=zeros(3*(n+ncm));                       %Creating the cubic part of the Global Mass Matrix
for i=1:3:3*n-2
    M1G(i:i+3,i:i+3)=M1G(i:i+3,i:i+3)+M1E;
end

%% M1G element 41 mapping: 1=97, 2=122,3=123,4=124

M1G(97,97)=M1G(97,97)+M1E(1,1);         %mapping column 1 of element 41
M1G(122,97)=M1G(122,97)+M1E(2,1);         %the nodes touch the rest of
M1G(123,97)=M1G(123,97)+M1E(3,1);         %the structure
M1G(124,97)=M1G(124,97)+M1E(4,1);         

j=2;
for i=122:124
    M1G(97,i)=M1G(97,i)+M1E(1,j);   %mapping the rest of element 41
    M1G(122,i)=M1G(122,i)+M1E(2,j);   
    M1G(123,i)=M1G(123,i)+M1E(3,j);   
    M1G(124,i)=M1G(124,i)+M1E(4,j);   
    j=j+1;
end

                                %mapping elements 42 through 47
for i=42:47
    
        M1G(3*i-2,3*i-2)=M1G(3*i-2,3*i-2)+M1E(1,1);
        M1G(3*i-1,3*i-2)=M1G(3*i-1,3*i-2)+M1E(2,1);
        M1G(3*i,3*i-2)=M1G(3*i,3*i-2)+M1E(3,1);
        M1G(3*i+1,3*i-2)=M1G(3*i+1,3*i-2)+M1E(4,1);
        
        M1G(3*i-2,3*i-1)=M1G(3*i-2,3*i-1)+M1E(1,2);
        M1G(3*i-1,3*i-1)=M1G(3*i-1,3*i-1)+M1E(2,2);
        M1G(3*i,3*i-1)=M1G(3*i,3*i-1)+M1E(3,2);
        M1G(3*i+1,3*i-1)=M1G(3*i+1,3*i-1)+M1E(4,2);
        
        M1G(3*i-2,3*i)=M1G(3*i-2,3*i)+M1E(1,3);
        M1G(3*i-1,3*i)=M1G(3*i-1,3*i)+M1E(2,3);
        M1G(3*i,3*i)=M1G(3*i,3*i)+M1E(3,3);
        M1G(3*i+1,3*i)=M1G(3*i+1,3*i)+M1E(4,3);
        
        M1G(3*i-2,3*i+1)=M1G(3*i-2,3*i+1)+M1E(1,4);
        M1G(3*i-1,3*i+1)=M1G(3*i-1,3*i+1)+M1E(2,4);
        M1G(3*i,3*i+1)=M1G(3*i,3*i+1)+M1E(3,4);
        M1G(3*i+1,3*i+1)=M1G(3*i+1,3*i+1)+M1E(4,4);
         
end

%%

% M1G element 48 mapping: 1=142, 2=143, 3=144, 4=25

j=1;
for i=142:144
    M1G(142,i)=M1G(142,i)+M1E(1,j);       %mapping the first 3 columns of
    M1G(143,i)=M1G(143,i)+M1E(2,j);       %element 48
    M1G(144,i)=M1G(144,i)+M1E(3,j);
    M1G(25,i)=M1G(25,i)+M1E(4,j);
    j=j+1;
end

M1G(142,25)=M1G(142,25)+M1E(1,4);           %mapping rest of element 48
M1G(143,25)=M1G(143,25)+M1E(2,4);           %these nodes touch the rest of 
M1G(144,25)=M1G(144,25)+M1E(3,4);           %the structure
M1G(25,25)=M1G(25,25)+M1E(4,4);    
%%
            %Creating Damping Matrix: D1G
D1G=(M1G/(rho*A))*mu;

Dst=zeros(5*(n+ncm));                       %Creating the Damping state Matrix
Dst(1:3*(n+ncm),1:3*(n+ncm))=D1G;

%%

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

 M2G=zeros(2*(n+ncm));                       %Creating the quadratic part of the Global Mass Matrix
for i=1:2:2*n-1
    M2G(i:i+2,i:i+2)=M2G(i:i+2,i:i+2)+M2E;
end

 %% M2G element 41 mapping: 1=65, 2=82,3=83
 
 M2G(65,65)=M2G(65,65)+M2E(1,1);         %mapping column 1 of element 41
 M2G(82,65)=M2G(82,65)+M2E(2,1);         %the nodes touch the rest of
 M2G(83,65)=M2G(83,65)+M2E(3,1);         %the structure
          

j=2;
for i=82:83
    M2G(65,i)=M2G(65,i)+M2E(1,j);   %mapping the rest of element 41
    M2G(82,i)=M2G(82,i)+M2E(2,j);   
    M2G(83,i)=M2G(83,i)+M2E(3,j);     
    j=j+1;
end
 
 %%  M2G element 48 mapping: 1=95, 2=96, 3=17
 
 j=1;
for i=95:96
    M2G(95,i)=M2G(95,i)+M2E(1,j);       %mapping the first 3 columns of
    M2G(96,i)=M2G(96,i)+M2E(2,j);       %element 48
    M2G(17,i)=M2G(17,i)+M2E(3,j);
    j=j+1;
end

M2G(95,17)=M2G(95,17)+M2E(1,3);           %mapping last column of element 48
M2G(96,17)=M2G(96,17)+M2E(2,3);           %these nodes touch the rest of 
M2G(17,17)=M2G(17,17)+M2E(3,3);           %the structure


for i=42:47                             %mapping elements 42 through 47
    
        M2G(2*i-1,2*i-1)=M2G(2*i-1,2*i-1)+M2E(1,1);
        M2G(2*i,2*i-1)=M2G(2*i,2*i-1)+M2E(2,1);
        M2G(2*i+1,2*i-1)=M2G(2*i+1,2*i-1)+M2E(3,1);
        
        M2G(2*i-1,2*i)=M2G(2*i-1,2*i)+M2E(1,2);
        M2G(2*i,2*i)=M2G(2*i,2*i)+M2E(2,2);
        M2G(2*i+1,2*i)=M2G(2*i+1,2*i)+M2E(3,2);
        
        M2G(2*i-1,2*i+1)=M2G(2*i-1,2*i+1)+M2E(1,3);
        M2G(2*i,2*i+1)=M2G(2*i,2*i+1)+M2E(2,3);
        M2G(2*i+1,2*i+1)=M2G(2*i+1,2*i+1)+M2E(3,3);
end

%%

Mst=zeros(5*(n+ncm));                       %Creating the Mass state Matrix
Mst(1:3*(n+ncm),1:3*(n+ncm))=M1G;
Mst(3*(n+ncm)+1:5*(n+ncm),3*(n+ncm)+1:5*(n+ncm))=M2G;
Mst=sparse(Mst);


%%

K1E=[];                        %Creating top left quadrant of element stiffness matrix
for i=1:4
    for j=1:4
     K1E(i,j)=(2/h)*G*A*ks*int(dphi(i)*dphi(j),-1,1);  
    end
end

K1G=zeros(3*(n+ncm));                       %Creating the cubic part of the Global stiffness Matrix
for i=1:3:3*n-2
    K1G(i:i+3,i:i+3)=K1G(i:i+3,i:i+3)+K1E;
end
%% 
                    
                    %K1G element 41 mapping: 1=97, 2=122,3=123,4=124
                    
K1G(97,97)=K1G(97,97)+K1E(1,1);         %mapping column 1 of element 41
K1G(122,97)=K1G(122,97)+K1E(2,1);         %the nodes touch the rest of
K1G(123,97)=K1G(123,97)+K1E(3,1);         %the structure
K1G(124,97)=K1G(124,97)+K1E(4,1);         

j=2;
for i=122:124
    K1G(97,i)=K1G(97,i)+K1E(1,j);   %mapping the rest of element 41
    K1G(122,i)=K1G(122,i)+K1E(2,j);   
    K1G(123,i)=K1G(123,i)+K1E(3,j);   
    K1G(124,i)=K1G(124,i)+K1E(4,j);   
    j=j+1;
end

                                %mapping elements 42 through 47
for i=42:47
    
        K1G(3*i-2,3*i-2)=K1G(3*i-2,3*i-2)+K1E(1,1);
        K1G(3*i-1,3*i-2)=K1G(3*i-1,3*i-2)+K1E(2,1);
        K1G(3*i,3*i-2)=K1G(3*i,3*i-2)+K1E(3,1);
        K1G(3*i+1,3*i-2)=K1G(3*i+1,3*i-2)+K1E(4,1);
        
        K1G(3*i-2,3*i-1)=K1G(3*i-2,3*i-1)+K1E(1,2);
        K1G(3*i-1,3*i-1)=K1G(3*i-1,3*i-1)+K1E(2,2);
        K1G(3*i,3*i-1)=K1G(3*i,3*i-1)+K1E(3,2);
        K1G(3*i+1,3*i-1)=K1G(3*i+1,3*i-1)+K1E(4,2);
        
        K1G(3*i-2,3*i)=K1G(3*i-2,3*i)+K1E(1,3);
        K1G(3*i-1,3*i)=K1G(3*i-1,3*i)+K1E(2,3);
        K1G(3*i,3*i)=K1G(3*i,3*i)+K1E(3,3);
        K1G(3*i+1,3*i)=K1G(3*i+1,3*i)+K1E(4,3);
        
        K1G(3*i-2,3*i+1)=K1G(3*i-2,3*i+1)+K1E(1,4);
        K1G(3*i-1,3*i+1)=K1G(3*i-1,3*i+1)+K1E(2,4);
        K1G(3*i,3*i+1)=K1G(3*i,3*i+1)+K1E(3,4);
        K1G(3*i+1,3*i+1)=K1G(3*i+1,3*i+1)+K1E(4,4);
         
end
%%

% K1G element 48 mapping: 1=142, 2=143, 3=144, 4=25

j=1;
for i=142:144
    K1G(142,i)=K1G(142,i)+K1E(1,j);       %mapping the first 3 columns of
    K1G(143,i)=K1G(143,i)+K1E(2,j);       %element 48
    K1G(144,i)=K1G(144,i)+K1E(3,j);
    K1G(25,i)=K1G(25,i)+K1E(4,j);
    j=j+1;
end

K1G(142,25)=K1G(142,25)+K1E(1,4);           %mapping rest of element 48
K1G(143,25)=K1G(143,25)+K1E(2,4);           %these nodes touch the rest of 
K1G(144,25)=K1G(144,25)+K1E(3,4);           %the structure
K1G(25,25)=K1G(25,25)+K1E(4,4);    

 %%             
  K2E=[];                        %Creating top right quadrant of element stiffness matrix
 for i=1:4
      for j=1:3
       K2E(i,j)=-G*A*ks*int(dphi(i)*g(j),-1,1);  
      end
 end 
  
 K2Expanded=[K2E(1:2,:);zeros(1,3);K2E(3:4,:)];
 K2Expanded=[K2Expanded(:,1),zeros(5,1),K2Expanded(:,2),zeros(5,1),K2Expanded(:,3)];
 
 %%
 %          K2G element 41 mapping:  phi1=129, phi2=162, phi3=164, phi4=165
                                    %  g1=129, g2=163, g3=165

 K2GExpanded=zeros(4*(n+ncm));
 for i=1:4:4*n-3
     K2GExpanded(i:i+4,i:i+4)=K2GExpanded(i:i+4,i:i+4)+K2Expanded;
 end

 K2GExpanded(129,129)=K2GExpanded(129,129)+K2E(1,1);
 K2GExpanded(162,129)=K2GExpanded(162,129)+K2E(2,1);
 K2GExpanded(164,129)=K2GExpanded(164,129)+K2E(3,1);
 K2GExpanded(165,129)=K2GExpanded(165,129)+K2E(4,1);
 
 j=2;
 for i=163:2:165
   K2GExpanded(129,i)=K2GExpanded(129,i)+K2E(1,j); 
   K2GExpanded(162,i)=K2GExpanded(162,i)+K2E(2,j);
   K2GExpanded(164,i)=K2GExpanded(164,i)+K2E(3,j);
   K2GExpanded(165,i)=K2GExpanded(165,i)+K2E(4,j);
   j=j+1;
 end

%  Mapping elements 42 to 47

for i=42:47
    
        K2GExpanded(4*i-3,4*i-3)=K2GExpanded(4*i-3,4*i-3)+K2Expanded(1,1);
        K2GExpanded(4*i-2,4*i-3)=K2GExpanded(4*i-2,4*i-3)+K2Expanded(2,1);
        K2GExpanded(4*i-1,4*i-3)=K2GExpanded(4*i-1,4*i-3)+K2Expanded(3,1);
        K2GExpanded(4*i,4*i-3)=K2GExpanded(4*i,4*i-3)+K2Expanded(4,1);
        K2GExpanded(4*i+1,4*i-3)=K2GExpanded(4*i+1,4*i-3)+K2Expanded(5,1);
        
        K2GExpanded(4*i-3,4*i-2)=K2GExpanded(4*i-3,4*i-2)+K2Expanded(1,2);
        K2GExpanded(4*i-2,4*i-2)=K2GExpanded(4*i-2,4*i-2)+K2Expanded(2,2);
        K2GExpanded(4*i-1,4*i-2)=K2GExpanded(4*i-1,4*i-2)+K2Expanded(3,2);
        K2GExpanded(4*i,4*i-2)=K2GExpanded(4*i,4*i-2)+K2Expanded(4,2);
        K2GExpanded(4*i+1,4*i-2)=K2GExpanded(4*i+1,4*i-2)+K2Expanded(5,2);
        
        K2GExpanded(4*i-3,4*i-1)=K2GExpanded(4*i-3,4*i-1)+K2Expanded(1,3);
        K2GExpanded(4*i-2,4*i-1)=K2GExpanded(4*i-2,4*i-1)+K2Expanded(2,3);
        K2GExpanded(4*i-1,4*i-1)=K2GExpanded(4*i-1,4*i-1)+K2Expanded(3,3);
        K2GExpanded(4*i,4*i-1)=K2GExpanded(4*i,4*i-1)+K2Expanded(4,3);
        K2GExpanded(4*i+1,4*i-1)=K2GExpanded(4*i+1,4*i-1)+K2Expanded(5,3);
        
        K2GExpanded(4*i-3,4*i)=K2GExpanded(4*i-3,4*i)+K2Expanded(1,4);
        K2GExpanded(4*i-2,4*i)=K2GExpanded(4*i-2,4*i)+K2Expanded(2,4);
        K2GExpanded(4*i-1,4*i)=K2GExpanded(4*i-1,4*i)+K2Expanded(3,4);
        K2GExpanded(4*i,4*i)=K2GExpanded(4*i,4*i)+K2Expanded(4,4);
        K2GExpanded(4*i+1,4*i)=K2GExpanded(4*i+1,4*i)+K2Expanded(5,4);
       
        K2GExpanded(4*i-3,4*i+1)=K2GExpanded(4*i-3,4*i+1)+K2Expanded(1,5);
        K2GExpanded(4*i-2,4*i+1)=K2GExpanded(4*i-2,4*i+1)+K2Expanded(2,5);
        K2GExpanded(4*i-1,4*i+1)=K2GExpanded(4*i-1,4*i+1)+K2Expanded(3,5);
        K2GExpanded(4*i,4*i+1)=K2GExpanded(4*i,4*i+1)+K2Expanded(4,5);
        K2GExpanded(4*i+1,4*i+1)=K2GExpanded(4*i+1,4*i+1)+K2Expanded(5,5);
        
         
end
 
 
% 
% %%K2G element 48 mapping:  phi1=189, phi2=190, phi3=192, phi4=33
%                      %  g1=189, g2=191, g3=33

j=1;
for i=189:2:191
    K2GExpanded(189,i)=K2GExpanded(189,i)+K2E(1,j);
    K2GExpanded(190,i)=K2GExpanded(190,i)+K2E(2,j);
    K2GExpanded(192,i)=K2GExpanded(192,i)+K2E(3,j);
    K2GExpanded(33,i)=K2GExpanded(33,i)+K2E(4,j);
    j=j+1;
end




K2GExpanded(189,33)=K2GExpanded(189,33)+K2E(1,3);
K2GExpanded(190,33)=K2GExpanded(190,33)+K2E(2,3);
K2GExpanded(192,33)=K2GExpanded(192,33)+K2E(3,3);
K2GExpanded(33,33)=K2GExpanded(33,33)+K2E(4,3);

 K2G=K2GExpanded(:,1:2:end);  %Deleting the even columns to get back to K2G
 remove=3:4:4*(n+ncm);              %Deleting every fourth row after row three to get back to K2G
 K2G(remove,:)=[];  
 
 %%
 K3G=K2G';                       %Creating bottom left quadrant of stiffness matrix
 
 %%
 
 
 
for i=1:3
    for j=1:3
     K4Ea(i,j)=E*I*(2/h)*int(dg(i)*dg(j),-1,1);
     K4Eb(i,j)=G*A*ks*(h/2)*int(g(i)*g(j),-1,1);
     K4E(i,j)=K4Ea(i,j)+K4Eb(i,j);
    end
end 


 
 K4G=zeros(2*(n+ncm));                       %Creating the bottom right quadrant of Global stiffness Matrix
for i=1:2:2*n-1
    K4G(i:i+2,i:i+2)=K4G(i:i+2,i:i+2)+K4E;
end
% 
% 
 %% K4G element 41 mapping: 1=65, 2=82,3=83
 
 K4G(65,65)=K4G(65,65)+K4E(1,1);         %mapping column 1 of element 41
 K4G(82,65)=K4G(82,65)+K4E(2,1);         %the nodes touch the rest of
 K4G(83,65)=K4G(83,65)+K4E(3,1);         %the structure
          

j=2;
for i=82:83
    K4G(65,i)=K4G(65,i)+K4E(1,j);   %mapping the rest of element 41
    K4G(82,i)=K4G(82,i)+K4E(2,j);   
    K4G(83,i)=K4G(83,i)+K4E(3,j);     
    j=j+1;
end
 
 %%  K4G element 48 mapping: 1=95, 2=96, 3=17
 
 j=1;
for i=95:96
    K4G(95,i)=K4G(95,i)+K4E(1,j);       %mapping the first 3 columns of
    K4G(96,i)=K4G(96,i)+K4E(2,j);       %element 48
    K4G(17,i)=K4G(17,i)+K4E(3,j);
    j=j+1;
end

K4G(95,17)=K4G(95,17)+K4E(1,3);           %mapping last column of element 48
K4G(96,17)=K4G(96,17)+K4E(2,3);           %these nodes touch the rest of 
K4G(17,17)=K4G(17,17)+K4E(3,3);           %the structure


for i=42:47                             %mapping elements 42 through 47
    
        K4G(2*i-1,2*i-1)=K4G(2*i-1,2*i-1)+K4E(1,1);
        K4G(2*i,2*i-1)=K4G(2*i,2*i-1)+K4E(2,1);
        K4G(2*i+1,2*i-1)=K4G(2*i+1,2*i-1)+K4E(3,1);
        
        K4G(2*i-1,2*i)=K4G(2*i-1,2*i)+K4E(1,2);
        K4G(2*i,2*i)=K4G(2*i,2*i)+K4E(2,2);
        K4G(2*i+1,2*i)=K4G(2*i+1,2*i)+K4E(3,2);
        
        K4G(2*i-1,2*i+1)=K4G(2*i-1,2*i+1)+K4E(1,3);
        K4G(2*i,2*i+1)=K4G(2*i,2*i+1)+K4E(2,3);
        K4G(2*i+1,2*i+1)=K4G(2*i+1,2*i+1)+K4E(3,3);
end

%%
 Kst=zeros(5*(n+ncm));               %Combining global stiffness matrix
 Kst(1:3*(n+ncm),1:3*(n+ncm))=K1G;
 Kst(1:3*(n+ncm),3*(n+ncm)+1:5*(n+ncm))=K2G;
 Kst(3*(n+ncm)+1:5*(n+ncm),1:3*(n+ncm))=K3G;
 Kst(3*(n+ncm)+1:5*(n+ncm),3*(n+ncm)+1:5*(n+ncm))=K4G;
 
 %%    
    
Kst(:,121)=[];                           %Applying boundary conditions to Stiffness
Kst(121,:)=[];
Kst(97,:)=[];
Kst(:,97)=[];
Kst(73,:)=[];
Kst(:,73)=[];
Kst(49,:)=[];
Kst(:,49)=[];
Kst(25,:)=[];
Kst(:,25)=[];
Kst(1,:)=[];
Kst(:,1)=[];

Mst(:,121)=[];                           %Applying boundary conditions to Mass
Mst(121,:)=[];
Mst(97,:)=[];
Mst(:,97)=[];
Mst(73,:)=[];
Mst(:,73)=[];
Mst(49,:)=[];
Mst(:,49)=[];
Mst(25,:)=[];
Mst(:,25)=[];
Mst(1,:)=[];
Mst(:,1)=[];

Dst(:,121)=[];                           %Applying boundary conditions to Damping
Dst(121,:)=[];
Dst(97,:)=[];
Dst(:,97)=[];
Dst(73,:)=[];
Dst(:,73)=[];
Dst(49,:)=[];
Dst(:,49)=[];
Dst(25,:)=[];
Dst(:,25)=[];
Dst(1,:)=[];
Dst(:,1)=[];


 Mst=sparse(Mst);
 Kst=sparse(Kst);
 
 q=-load*sin(2*pi*f*t);                       %Creating Force Element Vector
 ForceE=[];
 for i=1:4                                 
    ForceE(i)=int(phi(i),-1,1)*(h/2);
 end
 ForceE=ForceE';                             
 
 Qst=zeros(5*(n+ncm),1);                           %Creating the full Force State vector
 for i=1:3:3*n-2
     Qst(i:i+3)=Qst(i:i+3)+ForceE;
 end
                            
% 
 Qst(1:48,:)=0;      %Loading on the top of the structure only
 Qst(74:121,:)=0;
                            
                            
%                       %Applying boundary conditions to load     
Qst(121,:)=[];
Qst(97,:)=[];
Qst(73,:)=[];
Qst(49,:)=[];
Qst(25,:)=[];
Qst(1,:)=[];
 
 
 
  
 ui=zeros(5*(n+ncm)-6,1);                    %setting initial conditions ui=initial position state vector
 dui=zeros(5*(n+ncm)-6,1);                   %dui=initial velocity state vector
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
       COMP=[Mst+(Kst*(dt^2)*(1/4))+(Dst*(dt/2))];
      ddu=inv(COMP)*((q*Qst)-(Kst*(ui+dui*dt+(1/4)*ddui*dt^2))-Dst*(dui+ddui*dt*.5));
      du=dui+(.5*ddui+.5*ddu)*dt;
      u=ui+dui*dt+.5*(.5*ddui+.5*ddu)*dt^2;
      
      ui=u;           %setting state solution vector equal to initial for next time step
       dui=du;
       ddui=ddu;
 
        U(:,:,i)=u;                %3D matrix that holds displacement vectors for each time step
   end
  
 
 fixednode=zeros(1,1,Ntimestep);     
 U1=U(1:23,:,:);                %cutting just displacement vector from full state vector
 U2=U(24:46,:,:);
 U3=U(47:69,:,:);
 U4=U(70:92,:,:);
 U5=U(93:115,:,:);
 U6=U(116:138,:,:);
  
                                                 %putting fixed nodes back in displacement vector
  ULC=[fixednode;U1;fixednode;U2;fixednode];     %Left Column
  URC=[fixednode;U4;fixednode;U5;fixednode];     %Right column
  UTB=[fixednode;U3;fixednode];                  %Top Beam
  UMB=[fixednode;U6;fixednode];                   %Middle Beam
  y=0:48;        
  yh=y*(1/3)*h;                   %vector for column plot
  x=0:24;
  xh=x*(1/3)*h;                   %vector for top and mid beam plot
  z=0:48;
  zh=y*(1/3)*h;
  
  figure
  for i=1:Ntimestep
 URC(:,:,i)=-flipud(URC(:,:,i));
  x3=ones(49,1);
  y3=ones(25,1);
  y3v=ones(49,1);
  z3=ones(25,1);
 
  axis manual
  
  plot3(xh+2,y3,4*z3+10^6*UTB(:,:,i),xh+2,y3,2*z3-10^6*UMB(:,:,i),2*x3-10^6*ULC(:,:,i),y3v,zh,4*x3-10^6*URC(:,:,i),y3v,zh)
  axis([0 5 0 2 0 5 ])
  view(0,0);
  %legend('top','mid')
  xlabel('x (meters) all deformations magnified by a factor 10^6')
  %ylabel('y (meters)')
  zlabel('z (meters)')
  title([num2str(load),'N, 2 Hz distributed load on top beam.  t=',num2str(i/1000),' seconds'])
%  
  drawnow
  pause(1/120)
  end
