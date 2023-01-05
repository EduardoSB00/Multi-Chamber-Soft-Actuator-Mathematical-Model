P_o = 101325; %presion inicial del sistema en Pa
Q = 0.5; %Flujo volumetrico del compresor
W = zeros(1,100)
t = 1:100;

V_o = (7/1000)*(5/1000)*(11.5/1000);
for i=1:100
W(i) = P_o*V_o*Q*t(i);
end

plot(t,W, 'LineWidth',2);
xlabel('Tiempo (s)');
ylabel('Energia por Presion y Deformacion Volumetrica (J)');
grid on;
title 'Energia por Presion con Respecto al Tiempo';
% axis([0 0.08 -0.06 0.02]);
pause(0.1);