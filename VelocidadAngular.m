%% Definicion de Parametros
P_o = 101325; %presion inicial del sistema en Pa
Q = 0.5; %Flujo volumetrico del compresor
n = 7; %numero de camaras neumaticas
t=1:100; %Variable temporal
Mu = 18.3; %%Modulo de Corte en GPa

b=5.75/1000; %5.75 Longitud desde la base hasta el punto critico de la pared lateral
L = 11.5/1000; %Altura total de la camara 11.5
d = 5/1000; %Grosor de la camara 
e = 7/1000; %Ancho de la camara 7
c = 0.5/1000; %Distancia entre las camara 
z=10/1000; %Longitud de las camaras 


V_o = L*d*e;


vel = zeros(1,100); %Energia

W = zeros(1,100); %Energia
W1 = zeros(1,100); %Energia
Wt = zeros(1,100); %Energia
Rho = 1.225;
T = zeros(1,100); 
lambda = zeros(1,100);  %Deformacion lateral
h = zeros(1,100);
A = zeros(1,100); 
B = zeros(1,100); 
C = zeros(1,100); 
D = zeros(1,100); 
E = zeros(1,100); 
p = zeros(1,100); 
q = zeros(1,100); 
r = zeros(1,100);
p2 = zeros(1,100); 
q2 = zeros(1,100); 
Q2 = zeros(1,100); 
phi2 = zeros(1,100); 
beta0 = zeros(1,100); 
alpha0 = zeros(1,100); 
x = zeros(1,100); 
x0 = zeros(1,100); 
y = zeros(1,100); 
theta = zeros(1,100); 
Ox = zeros(1,100); 
Oy = zeros(1,100); 
thetaT = zeros(1,100);
thetaT(1) = 135;
lambda0 = zeros(1,100); 
h0 = zeros(1,100); 
A_t = (10^-6); 

Wk= zeros(1,100);

%% Computacion de Energia por Presion
for i=1:100
%%W(i) = P_o*V_o(i)*log10((V_o(i)+Q*t)/V_o(i) );%Dependiendo de la energia producida cambiara el angulo del actuador.
%%W(i)=k;
%%lambda(i)=1.0180;
W(i) = P_o*V_o*Q*t(i);
W1(i) = (P_o*pi*0.004*((lambda(i)*0.0005)^2))/n;%%cambiar ese valor de radio 0.00008
Wk(i)= P_o*V_o*Q*t(i);
end
%% Computacion deformacion lateral
for i=1:100
T(i) = sqrt(-1-(((-((2*W(i)/Mu) + 3))^3)/27));
lambda0(i)= (-2*nthroot(sqrt(((-1)^2)+(0^2)),3)*cos(180+(atan(0/-1))/3))-1;
lambda(i)= (-2*nthroot(sqrt(((-1)^2)+(T(i)^2)),3)*cos(180+(atan(T(i)/-1))/3))-lambda0(i); %%Lambda0 es para mapear correctamente los datos en la parte lineal del estiramiento.
%%h(i) = sqrt((2*(((2*lambda(i)*L(i))/(2*pi))^2))-((L(i)/2)^2)); 
%%prueba de aproximacion fallida por numeros imaginarios de gran magnitud.
h0(i) = sqrt(((((2*L)/2)^2)/2)-((L/2)^2));
h(i) = sqrt(((((2*lambda(i)*L)/2)^2)/2)-((L/2)^2))-h0(i);
end

for i=1:100
A(i) = (((b^2)-4*((h(i))^2))^2)+ (16*(b^2)*(h(i)^2));
B(i) = 32*c*b*(h(i))^2;
C(i) = 16*((c*h(i))^2)-16*(b^2)*(h(i)^2)+8*(h(i)^2)*((b^2)-4*(h(i)^2));
D(i) = -32*c*b*(h(i))^2;
E(i) = 16*(h(i)^4)-16*(c*h(i))^2;

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

for i=2:100
   thetaT(i) = thetaT(i-1) + theta(i); 
vel(i) = theta(i)/t(i);
end

plot(t,vel,'b', 'LineWidth',2);
xlabel('Tiempo (s)');
ylabel('Velocidad Angular (grad/s)');
grid on;
title 'Velocidad Angular Promedio con Respecto al Tiempo';
% axis([0 0.08 -0.06 0.02]);
pause(0.1);