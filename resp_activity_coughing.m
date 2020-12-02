%This script determines the expected rate of shedding SARS-CoV-2 via  
%droplets and aerosols while coughing. 

clear all
close all
clc

%Define values
p = 997e3; %Density of water at room temperature {g/m3}
Vml_Vg = 1; %Conversion factor (volume) from g to ml for water at room
            %temperature
v_SARS2_1daypre = 5.48; %SARS-CoV-2 mean viral load, -1 day {log10 copies/ml}
uCI_1daypre = 6.70; %upper 95% CI for v_SARS2_1daypre {log10 copies/ml}
lCI_1daypre = 4.25; %lower 95% CI for v_SARS2_1daypre {log10 copies/ml}
v_SARS2_1daypre_top10 = 10.01992; %SARS-CoV-2 viral load (98th cp), -1 day {log10 copies/ml}
uCI_1daypre_top10 = 12.51300;
lCI_1daypre_top10 = 8.02355;

v_SARS2_1daypost = 7.08; %SARS-CoV-2 mean viral load, -1 day {log10 copies/ml}
uCI_1daypost = 7.63;
lCI_1daypost = 6.54;
v_SARS2_1daypost_top10 = 11.35414; %SARS-CoV-2 viral load (98th cp), -1 day {log10 copies/ml}
uCI_1daypost_top10 = 12.35780;
lCI_1daypost_top10 = 10.43199;

gamma = 0.1e-2; %viability proportion {%}

%particle profile
d = [0.1653449
0.178063421
0.191760265
0.206510686
0.222395726
0.239502662
0.257925483
0.277765409
0.299131445
0.322140981
0.346920437
0.373605958
0.402344159
0.433292936
0.466622329
0.502515457
0.541169527
0.582796911
0.627626322
0.675904064
0.727895384
0.783885938
0.848663711
0.913943904
0.984245524
1.059954825
1.141487774
1.229292331
1.323850917
1.425683057
1.535348244
1.653449003
1.780634209
1.917602648
2.065106857
2.22395726
2.39502662
2.579254832
2.777654091
2.99131445
3.221409811
3.469204372
3.755888091
4.044795338
4.355925665
4.690988495
5.051824745
5.440416936
5.858900087
6.243128858
6.759040636
7.278953844
7.838859376
8.397266396
8.995451959
9.687392507
10.37748117
11.175729
12.03537898
13.02994339
14.03222167
15.0318175
16.18808225
25.79254832
40.44795338
64.10568596
102.1398939
164.4719926
259.2943797
404.4795338
644.4591621]*1e-6; %hydrated particle diameters for given respiratory activity {m}
V = (p*pi*Vml_Vg/6)* d.^3; %particle volumes {ml}
R = [5.33E-05
1.07E-04
1.77E-04
3.30E-04
5.62E-04
9.16E-04
1.56E-03
2.51E-03
3.80E-03
6.11E-03
9.38E-03
1.40E-02
2.12E-02
3.16E-02
4.25E-02
5.81E-02
7.92E-02
1.08E-01
1.39E-01
1.76E-01
2.20E-01
2.96E-01
3.70E-01
4.76E-01
6.21E-01
8.11E-01
1.04E+00
1.30E+00
1.51E+00
1.65E+00
1.70E+00
1.63E+00
1.45E+00
1.21E+00
9.99E-01
8.00E-01
6.31E-01
5.05E-01
4.10E-01
3.44E-01
2.75E-01
2.27E-01
1.82E-01
1.41E-01
1.10E-01
8.28E-02
6.07E-02
4.45E-02
3.12E-02
2.32E-02
1.51E-02
9.96E-03
6.48E-03
4.28E-03
2.83E-03
1.73E-03
1.11E-03
6.81E-04
4.37E-04
3.02E-04
2.60E-04
2.72E-04
3.69E-04
3.07E-03
1.54E-02
4.42E-02
7.21E-02
6.60E-02
3.39E-02
1.02E-02
1.72E-03]; %rate of particle being expelled, corresponding to particle diameters {particles/interval}


%Initialize matrices for Poisson lamda, partitioning likelihood, variance 
%(via delta method) and expelling rate at each size
n = [0:1:140]'; %up to 140 viable virions

lambda_SARS2_1daypre = zeros(length(d),1);
uCI_1daypre_lambda = zeros(length(d),1);
lCI_1daypre_lambda = zeros(length(d),1);
P_SARS2_1daypre = zeros(length(d),length(n));
Var_SARS2_1daypre = zeros(length(d),length(n));

lambda_SARS2_1daypre_top10 = zeros(length(d),1);
uCI_1daypre_top10_lambda = zeros(length(d),1);
lCI_1daypre_top10_lambda = zeros(length(d),1);
P_SARS2_1daypre_top10 = zeros(length(d),length(n));
Var_SARS2_1daypre_top10 = zeros(length(d),length(n));

lambda_SARS2_1daypost = zeros(length(d),1);
uCI_1daypost_lambda = zeros(length(d),1);
lCI_1daypost_lambda = zeros(length(d),1);
P_SARS2_1daypost = zeros(length(d),length(n));
Var_SARS2_1daypost = zeros(length(d),length(n));

lambda_SARS2_1daypost_top10 = zeros(length(d),1);
uCI_1daypost_top10_lambda = zeros(length(d),1);
lCI_1daypost_top10_lambda = zeros(length(d),1);
P_SARS2_1daypost_top10 = zeros(length(d),length(n));
Var_SARS2_1daypost_top10 = zeros(length(d),length(n));

Expel = zeros(length(d),4); %Rate of viable SARS-CoV-2 virions being expelled
uCI = zeros(length(d),4); %upper 95% CI for each expel
lCI = zeros(length(d),4); %lower 95% CI for each expel


%Determine rates of expelling SARS-CoV-2 for the given respiratory activity
for i = 1 : length(d)
    lambda_SARS2_1daypre(i) = (10^v_SARS2_1daypre)* V(i)*gamma;
    uCI_1daypre_lambda(i) = (10^uCI_1daypre)* V(i)*gamma;
    lCI_1daypre_lambda(i) = (10^lCI_1daypre)* V(i)*gamma;
    lambda_SARS2_1daypre_top10(i) = (10^v_SARS2_1daypre_top10)* V(i)*gamma;
    uCI_1daypre_top10_lambda(i) = (10^uCI_1daypre_top10)* V(i)*gamma;
    lCI_1daypre_top10_lambda(i) = (10^lCI_1daypre_top10)* V(i)*gamma;
    lambda_SARS2_1daypost(i) = (10^v_SARS2_1daypost)* V(i)*gamma;
    uCI_1daypost_lambda(i) = (10^uCI_1daypost)* V(i)*gamma;
    lCI_1daypost_lambda(i) = (10^lCI_1daypost)* V(i)*gamma;
    lambda_SARS2_1daypost_top10(i) = (10^v_SARS2_1daypost_top10)* V(i)*gamma;
    uCI_1daypost_top10_lambda(i) = (10^uCI_1daypost_top10)* V(i)*gamma;
    lCI_1daypost_top10_lambda(i) = (10^lCI_1daypost_top10)* V(i)*gamma;
    
        for k = 1 : length(n)
        P_SARS2_1daypre(i,k) = (lambda_SARS2_1daypre(i)^n(k))* ...
            exp(-lambda_SARS2_1daypre(i))/ factorial(n(k));
        Var_SARS2_1daypre(i,k) = lambda_SARS2_1daypre(i)* ...
            ( exp(-lambda_SARS2_1daypre(i))* lambda_SARS2_1daypre(i)^...
            (n(k)-1)* (n(k) - lambda_SARS2_1daypre(i))/ factorial(n(k)) )^2;
        
        P_SARS2_1daypre_top10(i,k) = (lambda_SARS2_1daypre_top10(i)^n(k))* ...
            exp(-lambda_SARS2_1daypre_top10(i))/ factorial(n(k));
        Var_SARS2_1daypre_top10(i,k) = lambda_SARS2_1daypre_top10(i)* ...
            ( exp(-lambda_SARS2_1daypre_top10(i))* lambda_SARS2_1daypre_top10(i)^...
            (n(k)-1)* (n(k) - lambda_SARS2_1daypre_top10(i))/ factorial(n(k)) )^2;
        
        P_SARS2_1daypost(i,k) = (lambda_SARS2_1daypost(i)^n(k))* ...
            exp(-lambda_SARS2_1daypost(i))/ factorial(n(k));
        Var_SARS2_1daypost(i,k) = lambda_SARS2_1daypost(i)* ...
            ( exp(-lambda_SARS2_1daypost(i))* lambda_SARS2_1daypost(i)^...
            (n(k)-1)* (n(k) - lambda_SARS2_1daypost(i))/ factorial(n(k)) )^2;
        
        P_SARS2_1daypost_top10(i,k) = (lambda_SARS2_1daypost_top10(i)^n(k))* ...
            exp(-lambda_SARS2_1daypost_top10(i))/ factorial(n(k));
        Var_SARS2_1daypost_top10(i,k) = lambda_SARS2_1daypost_top10(i)* ...
            ( exp(-lambda_SARS2_1daypost_top10(i))* lambda_SARS2_1daypost_top10(i)^...
            (n(k)-1)* (n(k) - lambda_SARS2_1daypost_top10(i))/ factorial(n(k)) )^2;
        
        end
end

Expel(:,1) = lambda_SARS2_1daypre.*R;
Expel(:,2) = lambda_SARS2_1daypre_top10.*R;
Expel(:,3) = lambda_SARS2_1daypost.*R;
Expel(:,4) = lambda_SARS2_1daypost_top10.*R;

%Since the above is equivalent to using the expectation value multiplied by
%the particle rate, variances can be calculated using lambda
uCI(:,1) = uCI_1daypre_lambda.*R - Expel(:,1);
uCI(:,2) = uCI_1daypre_top10_lambda.*R - Expel(:,2);
uCI(:,3) = uCI_1daypost_lambda.*R - Expel(:,3);
uCI(:,4) = uCI_1daypost_top10_lambda.*R - Expel(:,4);

lCI(:,1) = Expel(:,1) - lCI_1daypre_lambda.*R;
lCI(:,2) = Expel(:,2) - lCI_1daypre_top10_lambda.*R;
lCI(:,3) = Expel(:,3) - lCI_1daypost_lambda.*R;
lCI(:,4) = Expel(:,4) - lCI_1daypost_top10_lambda.*R;

%{
Var(:,1) = lamda_SARS2_1daypre;
Var(:,2) = lamda_SARS2_1daypre_top10;
Var(:,3) = lamda_SARS2_1daypost;
Var(:,4) = lamda_SARS2_1daypost_top10;

for i = 1:length(d)
    Var(i,1) = sum(P_SARS2_1daypre(i,:)*n*...
        Var_SARS2_1daypre(i,:))/sum(P_SARS2_1daypre(i,:));
    Var(i,2) = sum(P_SARS2_1daypre_top10(i,:)*n*...
        Var_SARS2_1daypre_top10(i,:))/sum(P_SARS2_1daypre_top10(i,:));
    Var(i,3) = sum(P_SARS2_1daypost(i,:)*n*...
        Var_SARS2_1daypost(i,:))/sum(P_SARS2_1daypost(i,:));
    Var(i,4) = P_SARS2_1daypost_top10(i,:).*n'*...
        Var_SARS2_1daypost_top10(i,:)'/sum(P_SARS2_1daypost_top10(i,:));
end


Var2 = zeros(length(d),4); %Variance

for i = 1:length(d)
    Var2(i,1) = sum(P_SARS2_1daypre(i,2:end).*...
        Var_SARS2_1daypre(i,2:end))/sum(P_SARS2_1daypre(i,2:end));
    Var2(i,2) = sum(P_SARS2_1daypre_top10(i,2:end).*...
        Var_SARS2_1daypre_top10(i,2:end))/sum(P_SARS2_1daypre_top10(i,2:end));
    Var2(i,3) = sum(P_SARS2_1daypost(i,2:end).*...
        Var_SARS2_1daypost(i,2:end))/sum(P_SARS2_1daypost(i,2:end));
    Var2(i,4) = sum(P_SARS2_1daypost_top10(i,2:end).*...
        Var_SARS2_1daypost_top10(i,2:end))/sum(P_SARS2_1daypost_top10(i,2:end));
end
%}

subplot(2,2,1)
plot(d,Expel(:,1))
title('SARS-CoV-2, 1 day pre-symptom onset, mean')
set(gca,'YScale','log')
hold on
hold off

subplot(2,2,2)
plot(d,Expel(:,2))
title('SARS-CoV-2, 1 day pre-symptom onset, top 10%')
set(gca,'YScale','log')

subplot(2,2,3)
plot(d,Expel(:,3))
title('SARS-CoV-2, 1 day post-symptom onset, mean')
set(gca,'YScale','log')

subplot(2,2,4)
plot(d,Expel(:,4))
title('SARS-CoV-2, 1 day post-symptom onset, top 10%')
set(gca,'YScale','log')