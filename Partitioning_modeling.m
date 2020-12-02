%This script determines the likelihoods of respiratory droplets and 
%aerosols containing SARS-CoV-2 or influenza A(H1N1)pdm09 upon 
%atomization from the respiratory tract. It uses the delta method to 
%determine variance + various estimates of respiratory viral load. 

clear all
close all
clc

%Define values
p = 997e3; %Density of water at room temperature {g/m3}
Vml_Vg = 1; %Conversion factor (volume) from g to ml for water at room
            %temperature
v_SARS2_1daypre = 5.48; %SARS-CoV-2 mean viral load, -1 day {log10 copies/ml}
v_SARS2_1daypre_98 = 10.01992; %SARS-CoV-2 viral load (95th percentile), -1 day {log10 copies/ml}

v_SARS2_1daypost = 7.08; %SARS-CoV-2 mean viral load, -1 day {log10 copies/ml}
v_SARS2_1daypost_98 = 11.35414; %SARS-CoV-2 viral load (95th percentile), -1 day {log10 copies/ml}

v_H1N1 = 6.81; %A(H1N1)pdm09 mean viral load {log10 copies/ml}
v_H1N1_98 = 9.49969; %A(H1N1)pdm09 viral load (95th percentile), top 10% {log10 copies/ml}

gamma = 0.1e-2; %viability proportion {%}

d = linspace(0.1,100,334)'*1e-6; %particle diameters {m}
V = (p*pi*Vml_Vg/6)* d.^3; %particle volumes {ml}


%Initialize matrices for Poisson lamda, partitioning likelihood and
%variance (via delta method)
n = [0:1:140]; %up to 140 viable virions

lamda_SARS2_1daypre = zeros(length(d),1);
P_SARS2_1daypre = zeros(length(d),length(n));
Var_SARS2_1daypre = zeros(length(d),length(n));

lamda_SARS2_1daypre_98 = zeros(length(d),1);
P_SARS2_1daypre_98 = zeros(length(d),length(n));
Var_SARS2_1daypre_98 = zeros(length(d),length(n));

lamda_SARS2_1daypost = zeros(length(d),1);
P_SARS2_1daypost = zeros(length(d),length(n));
Var_SARS2_1daypost = zeros(length(d),length(n));

lamda_SARS2_1daypost_98 = zeros(length(d),1);
P_SARS2_1daypost_98 = zeros(length(d),length(n));
Var_SARS2_1daypost_98 = zeros(length(d),length(n));

lamda_H1N1 = zeros(length(d),1);
P_H1N1 = zeros(length(d),length(n));
Var_H1N1 = zeros(length(d),length(n));

lamda_H1N1_98 = zeros(length(d),1);
P_H1N1_98 = zeros(length(d),length(n));
Var_H1N1_98 = zeros(length(d),length(n));


%Determine likelihoods for virus partitioning based on Poisson statistics
for i = 1 : length(d)
    lamda_SARS2_1daypre(i) = (10^v_SARS2_1daypre)* V(i)*gamma;
    lamda_SARS2_1daypre_98(i) = (10^v_SARS2_1daypre_98)* V(i)*gamma;
    
    lamda_SARS2_1daypost(i) = (10^v_SARS2_1daypost)* V(i)*gamma;
    lamda_SARS2_1daypost_98(i) = (10^v_SARS2_1daypost_98)* V(i)*gamma;
    
    lamda_H1N1(i) = (10^v_H1N1)* V(i)*gamma;
    lamda_H1N1_98(i) = (10^v_H1N1_98)* V(i)*gamma;
    
        for k = 1 : length(n)
        P_SARS2_1daypre(i,k) = (lamda_SARS2_1daypre(i)^n(k))* ...
            exp(-lamda_SARS2_1daypre(i))/ factorial(n(k));
        Var_SARS2_1daypre(i,k) = lamda_SARS2_1daypre(i)* ...
            ( exp(-lamda_SARS2_1daypre(i))* lamda_SARS2_1daypre(i)^...
            (n(k)-1)* (n(k) - lamda_SARS2_1daypre(i))/ factorial(n(k)) )^2;
        
        P_SARS2_1daypre_98(i,k) = (lamda_SARS2_1daypre_98(i)^n(k))* ...
            exp(-lamda_SARS2_1daypre_98(i))/ factorial(n(k));
        Var_SARS2_1daypre_98(i,k) = lamda_SARS2_1daypre_98(i)* ...
            ( exp(-lamda_SARS2_1daypre_98(i))* lamda_SARS2_1daypre_98(i)^...
            (n(k)-1)* (n(k) - lamda_SARS2_1daypre_98(i))/ factorial(n(k)) )^2;
        
        P_SARS2_1daypost(i,k) = (lamda_SARS2_1daypost(i)^n(k))* ...
            exp(-lamda_SARS2_1daypost(i))/ factorial(n(k));
        Var_SARS2_1daypost(i,k) = lamda_SARS2_1daypost(i)* ...
            ( exp(-lamda_SARS2_1daypost(i))* lamda_SARS2_1daypost(i)^...
            (n(k)-1)* (n(k) - lamda_SARS2_1daypost(i))/ factorial(n(k)) )^2;
        
        P_SARS2_1daypost_98(i,k) = (lamda_SARS2_1daypost_98(i)^n(k))* ...
            exp(-lamda_SARS2_1daypost_98(i))/ factorial(n(k));
        Var_SARS2_1daypost_98(i,k) = lamda_SARS2_1daypost_98(i)* ...
            ( exp(-lamda_SARS2_1daypost_98(i))* lamda_SARS2_1daypost_98(i)^...
            (n(k)-1)* (n(k) - lamda_SARS2_1daypost_98(i))/ factorial(n(k)) )^2;
        
        P_H1N1(i,k) = (lamda_H1N1(i)^n(k))* exp(-lamda_H1N1(i)) /...
            factorial(n(k));
        Var_H1N1(i,k) = lamda_H1N1(i)* ...
            ( exp(-lamda_H1N1(i))* lamda_H1N1(i)^...
            (n(k)-1)* (n(k) - lamda_H1N1(i))/ factorial(n(k)) )^2;
        
        P_H1N1_98(i,k) = (lamda_H1N1_98(i)^n(k))* ...
            exp(-lamda_H1N1_98(i))/ factorial(n(k));
        Var_H1N1_98(i,k) = lamda_H1N1_98(i)* ...
            ( exp(-lamda_H1N1_98(i))* lamda_H1N1_98(i)^...
            (n(k)-1)* (n(k) - lamda_H1N1_98(i))/ factorial(n(k)) )^2;
        end
        
end

subplot(2,3,1)
plot(d*1e6, P_SARS2_1daypre(:,1))
hold on
plot(d*1e6, P_SARS2_1daypre(:,2))
plot(d*1e6, P_SARS2_1daypre(:,3))
plot(d*1e6, P_SARS2_1daypre(:,4))
plot(d*1e6, P_SARS2_1daypre(:,5))
hold off
title('SARS-CoV-2, -1 day (mean)')
xlabel('Hydrated diameter (um)')
ylabel('Likelihood of containining viable virion')
legend('empty','1 copy', '...')

subplot(2,3,2)
plot(d*1e6, P_SARS2_1daypre_98(:,1))
hold on
plot(d*1e6, P_SARS2_1daypre_98(:,2))
plot(d*1e6, P_SARS2_1daypre_98(:,3))
plot(d*1e6, P_SARS2_1daypre_98(:,4))
plot(d*1e6, P_SARS2_1daypre_98(:,5))
plot(d*1e6, P_SARS2_1daypre_98(:,6))
plot(d*1e6, P_SARS2_1daypre_98(:,7))
plot(d*1e6, P_SARS2_1daypre_98(:,10))
plot(d*1e6, P_SARS2_1daypre_98(:,15))
plot(d*1e6, P_SARS2_1daypre_98(:,20))
plot(d*1e6, P_SARS2_1daypre_98(:,30))
hold off
title('SARS-CoV-2, -1 day (98th percentile)')
xlabel('Hydrated diameter (um)')
ylabel('Likelihood of containining viable virion')
legend('empty','1 copy', '...')

subplot(2,3,3)
plot(d*1e6, P_SARS2_1daypost(:,1))
hold on
plot(d*1e6, P_SARS2_1daypost(:,2))
plot(d*1e6, P_SARS2_1daypost(:,3))
plot(d*1e6, P_SARS2_1daypost(:,4))
plot(d*1e6, P_SARS2_1daypost(:,5))
hold off
title('SARS-CoV-2, +1 day (mean)')
xlabel('Hydrated diameter (um)')
ylabel('Likelihood of containining viable virion')
legend('empty','1 copy', '...')

subplot(2,3,4)
plot(d*1e6, P_SARS2_1daypost_98(:,1))
hold on
plot(d*1e6, P_SARS2_1daypost_98(:,2))
plot(d*1e6, P_SARS2_1daypost_98(:,3))
plot(d*1e6, P_SARS2_1daypost_98(:,4))
plot(d*1e6, P_SARS2_1daypost_98(:,5))
plot(d*1e6, P_SARS2_1daypost_98(:,6))
plot(d*1e6, P_SARS2_1daypost_98(:,8))
plot(d*1e6, P_SARS2_1daypost_98(:,11))
plot(d*1e6, P_SARS2_1daypost_98(:,16))
plot(d*1e6, P_SARS2_1daypost_98(:,21))
plot(d*1e6, P_SARS2_1daypost_98(:,31))
plot(d*1e6, P_SARS2_1daypost_98(:,41))
plot(d*1e6, P_SARS2_1daypost_98(:,51))
plot(d*1e6, P_SARS2_1daypost_98(:,66))
plot(d*1e6, P_SARS2_1daypost_98(:,81))
plot(d*1e6, P_SARS2_1daypost_98(:,101))
plot(d*1e6, P_SARS2_1daypost_98(:,121))
plot(d*1e6, P_SARS2_1daypost_98(:,141))
hold off
title('SARS-CoV-2, +1 day (98th percentile)')
xlabel('Hydrated diameter (um)')
ylabel('Likelihood of containining viable virion')
legend('empty','1 copy', '...')

subplot(2,3,5)
plot(d*1e6, P_H1N1(:,1))
hold on
plot(d*1e6, P_H1N1(:,2))
plot(d*1e6, P_H1N1(:,3))
plot(d*1e6, P_H1N1(:,4))
plot(d*1e6, P_H1N1(:,5))
hold off
title('A(H1N1)pdm09, overall (mean)')
xlabel('Hydrated diameter (um)')
ylabel('Likelihood of containining viable virion')
legend('empty','1 copy', '...')

subplot(2,3,6)
plot(d*1e6, P_H1N1_98(:,1))
hold on
plot(d*1e6, P_H1N1_98(:,2))
plot(d*1e6, P_H1N1_98(:,3))
plot(d*1e6, P_H1N1_98(:,4))
plot(d*1e6, P_H1N1(:,5))
plot(d*1e6, P_H1N1(:,6))
plot(d*1e6, P_H1N1(:,7))
plot(d*1e6, P_H1N1(:,10))
hold off
title('A(H1N1)pdm09, overall (98th percentile)')
xlabel('Hydrated diameter (um)')
ylabel('Likelihood of containining viable virion')
legend('empty','1 copy', '...')