%% Definicion de Parametros
P_o = 101325; %presion inicial del sistema
b_1 = 7.25; %distancia desde la base hasta el punto critico en la pared lateral 1
L_1 = 14.5; %Altura total de la camara 1
d_1 = 11; %Grosor de la camara 1
e_1 = 7; %Ancho de la camara 1
V_o1 = L_1*d_1*e_1; %Volumen inicial de la camara 1
Q = 3; %Flujo volumetrico del compresor

%% Computacion de Energia por Presion
W = P_o*V_o*log((V_o+Qt)/V_o );