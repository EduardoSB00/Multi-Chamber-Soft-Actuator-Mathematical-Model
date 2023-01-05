%% Definicion de Parametros
P_o = 101325; %presion inicial del sistema en Pa
Qf = 2*0.0018877898; %Flujo volumetrico de referencia del compresor en m^3/s
n = 99; %numero de camaras neumaticas
t = linspace(0,10); %Variable temporal
Mu = 18000000000; %Modulo de Corte en Pa

Ang = 0; %Angulo entre dos camaras
 
b = zeros(1,n+1) + 5.75/1000; %5.75 Longitud desde la base hasta el punto critico de la pared lateral
L = zeros(1,n+1) + 11.5/1000; %Altura total de la camara 11.5
d = zeros(1,n+1) + 5/1000; %Grosor de la camara 
e = zeros(1,n+1) + 7/1000; %Ancho de la camara 7
c = zeros(1,n+1) + 0.5/1000; %Distancia entre las camara 
Z = zeros(1,n+1) + 10/1000; %Longitud de las camaras 

V_o = zeros(1,n+1)%Volumen inicial de las camaras 

for i=1:n+1
V_o(i) = L(i)*d(i)*e(i);
end

W = zeros(1,n+1); %Energia por Parte de la Presion
T = zeros(1,n+1); %Componente del Metodo de Cardano
lambda = zeros(1,n+1);  %Deformacion Lateral Unitaria
h = zeros(1,n+1); %Deformacion Radial Lateral

%%Parametros de Resolucion para el  Metodo de Ferrari
A = zeros(1,n+1); 
B = zeros(1,n+1); 
C = zeros(1,n+1); 
D = zeros(1,n+1); 
E = zeros(1,n+1); 
p = zeros(1,n+1); 
q = zeros(1,n+1); 
r = zeros(1,n+1);
p2 = zeros(1,n+1); 
q2 = zeros(1,n+1); 
Q2 = zeros(1,n+1); 
phi2 = zeros(1,n+1); 
beta0 = zeros(1,n+1); 
alpha0 = zeros(1,n+1); 
x = zeros(1,n+1); 
x0 = zeros(1,n+1); 
y = zeros(1,n+1); 

theta = zeros(1,n+1); %Angulo Relativo de Desplazamiento por Camara
Ox = zeros(1,n+1);  %Componente X del Vector de Direccion
Oy = zeros(1,n+1); %Componente Y del Vector de Direccion
thetaT = zeros(1,n+1); %Angulo Global de Desplazamiento por Camara
thetaT(1) = 135; %Angulo Inicial
lambda0 = zeros(1,n+1); %Error Inicial de Lambda por Aproximacion
h0 = zeros(1,n+1); %Error Inicial de h por Aproximacion

%for k= 0:100

%% Computacion de Energia por Presion
for i=1:n+1    
W(i) = (P_o*(((V_o(i)+(Qf)*t(i))^2)-(V_o(i)^2)))/(2*V_o(i));

end

%% Computacion deformacion lateral
for i=1:n+1
T(i) = sqrt(-1-(((-((2*W(i)/Mu) + 3))^3)/27));
lambda0(i)= (-2*nthroot(sqrt(((-1)^2)+(0^2)),3)*cos(180+(atan(0/-1))/3))-1;
lambda(i)= (-2*nthroot(sqrt(((-1)^2)+(T(i)^2)),3)*cos(180+(atan(T(i)/-1))/3))-lambda0(i); 
h0(i) = sqrt(((((2*L(i))/2)^2)/2)-((L(i)/2)^2));
h(i) = sqrt(((((2*lambda(i)*L(i))/2)^2)/2)-((L(i)/2)^2))-h0(i);
end

%% Computacion angulo de desplazamiento
for i=3:n+1
    
A(i) = (((b(i)^2)-4*((h(i))^2))^2)+ (16*(b(i)^2)*(h(i)^2));
B(i) = 32*c(i)*b(i)*(h(i))^2;
C(i) = 16*((c(i)*h(i))^2)-16*(b(i)^2)*(h(i)^2)+8*(h(i)^2)*((b(i)^2)-4*(h(i)^2));
D(i) = -32*c(i)*b(i)*(h(i))^2;
E(i) = 16*(h(i)^4)-16*(c(i)*h(i))^2;

p(i) = ((B(i)^2)/(8*A(i)^2))+(C(i)/A(i))-((B(i)^2)/(2*(A(i)^2)));
q(i) = ((B(i)^3)/(8*(A(i)^3)))-((C(i)*B(i))/(2*(A(i)^2)))+(D(i)/A(i));
r(i) = ((B(i)^4)/(256*A(i)^4))-(B(i)^4/(64*A(i)^4))+((C(i)*B(i)^2)/(16*A(i)^3))-((D(i)*B(i))/(4*A(i)^2))+(E(i)/A(i));

p2(i) = ((-8*p(i)^2/9)-(16*p(i)^2/9)+(16*p(i)^2/3)-2*p(i)^2+8*r(i))/-8;
q2(i) =(((2*p(i)^3)/3)-((8*p(i)*r(i))/3)+q(i)^2-((8*p(i)^3)/9)+((8*p(i)^3)/27))/-8;

Q2(i) =sqrt(-((q2(i)^2)/4)-((p2(i)^3)/27));
phi2(i) =atan(Q2(i)/-(q2(i)/2));

if (q2(i)^2/4)+(p2(i)^3/27) >= 0
    beta0(i) = nthroot(-(q2(i)/2)+sqrt((q2(i)^2/4)+(p2(i)^3/27)),3)+nthroot(-(q2(i)/2)-sqrt((q2(i)^2/4)+(p2(i)^3/27)),3);
else
    beta0(i) = 2*nthroot(sqrt((-q2(i)/2)^2+Q2(i)^2),3)*cos((180+phi2(i))/3);
end

alpha0(i) = beta0(i)-(p(i)/3);
x0(i)= q(i)/4*alpha0(i);

x(i) = (sqrt(abs(2*alpha0(i)))+sqrt(abs(2*alpha0(i)-4*(x0(i)*sqrt(abs(2*alpha0(i)))+(p(i)/2)+alpha0(i)))))/2;

y(i) = x(i) - B(i)/(4*A(i));

theta(i) = (180*asin(y(i)))/pi; 
end

for i=2:n+1
   thetaT(i) = thetaT(i-1) + theta(i); 
   Ox(i) = Ox(i-1) + Z(i-1)*cos((thetaT(i-1) + thetaT(i))*pi/180);
   Oy(i) = Oy(i-1) + Z(i-1)*sin((thetaT(i-1) + thetaT(i))*pi/180);
end

Ang = theta(3)*2;

%hold on;
t=t+0.1;
plot(t,h, 'LineWidth',2);
xlabel('Time (s)');
ylabel('Lateral Deformation (mm)');
grid on;
title 'Change in Deformation Dependent on Time';
%%axis([0 0.09 -0.07 0.02]);
pause(0.05);


%end