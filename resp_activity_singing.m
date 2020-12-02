%This script determines the expected rate of shedding SARS-CoV-2 via  
%droplets and aerosols while singing.

%Currently, data for talking (normally) are input for the particle profile

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
d = [0.276687125
0.298633223
0.322320029
0.347885609
0.375478985
0.405260995
0.437405236
0.472099075
0.509544737
0.549960492
0.593581918
0.640663282
0.691479016
0.746325322
0.805521892
0.86941378
0.938373405
1.012802727
1.07688963
1.156428991
1.252954263
1.351503671
1.457804348
1.565478174
1.673634257
1.805271676
1.938609532
2.081795759
2.255559911
2.411392505
2.577991294
2.768402431
2.986147347
3.235396275
3.474363523
3.714401076
3.988747749
4.264323186
4.599727708
4.939464846
5.327971715
5.721497153
6.116785143
6.568572556
7.116840768
7.642492561
8.206969139
8.813138109
9.464078885
10.16309832
10.91374753
11.77215345
12.75475621
13.75796416
14.77413113
15.79484938
16.88608723
17.97249338
19.38609532
20.63334536
22.05886624
23.16646573
24.22156179
26.97361948
42.64167484
66.72819729
107.1110444
171.0602771
270.4233565
425.3331983
679.2727599]*1e-6; %hydrated particle diameters for given respiratory activity {m}
V = (p*pi*Vml_Vg/6)* d.^3; %particle volumes {ml}
R = [3.7997E-04
7.0574E-04
1.3294E-03
2.4007E-03
4.3966E-03
7.1945E-03
1.1940E-02
2.0670E-02
3.1525E-02
4.9455E-02
8.0927E-02
1.1833E-01
2.0671E-01
3.3011E-01
5.5308E-01
9.3629E-01
1.4896E+00
2.2402E+00
2.1116E+01
2.6796E+01
3.2778E+01
3.8651E+01
4.3935E+01
4.7267E+01
5.0853E+01
5.0837E+01
4.9899E+01
4.7213E+01
4.3858E+01
4.0003E+01
3.7851E+01
3.5814E+01
3.3269E+01
3.2060E+01
3.1469E+01
2.9775E+01
2.7661E+01
2.6172E+01
2.3871E+01
2.2176E+01
2.0982E+01
1.9852E+01
1.7778E+01
1.5631E+01
1.2310E+01
1.0244E+01
8.3692E+00
7.0933E+00
6.1232E+00
5.5849E+00
5.1883E+00
4.8197E+00
4.4772E+00
3.8648E+00
3.2756E+00
2.6763E+00
2.2271E+00
1.7540E+00
1.2603E+00
9.2239E-01
6.5074E-01
4.5077E-01
3.1804E-01
9.0299E-04
6.3740E-03
2.3582E-02
4.9685E-02
5.4790E-02
3.3455E-02
1.1153E-02
1.9735E-03]; %rate of particle being expelled, corresponding to particle diameters {particles/interval}


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

        P_SARS2_1daypre_top10(i,k) = (lambda_SARS2_1daypre_top10(i)^n(k))* ...
            exp(-lambda_SARS2_1daypre_top10(i))/ factorial(n(k));
        
        P_SARS2_1daypost(i,k) = (lambda_SARS2_1daypost(i)^n(k))* ...
            exp(-lambda_SARS2_1daypost(i))/ factorial(n(k));
        
        P_SARS2_1daypost_top10(i,k) = (lambda_SARS2_1daypost_top10(i)^n(k))* ...
            exp(-lambda_SARS2_1daypost_top10(i))/ factorial(n(k));
        
        end
end

Expel(:,1) = lambda_SARS2_1daypre.*R;
Expel(:,2) = lambda_SARS2_1daypre_top10.*R;
Expel(:,3) = lambda_SARS2_1daypost.*R;
Expel(:,4) = lambda_SARS2_1daypost_top10.*R;

uCI(:,1) = uCI_1daypre_lambda.*R - Expel(:,1);
uCI(:,2) = uCI_1daypre_top10_lambda.*R - Expel(:,2);
uCI(:,3) = uCI_1daypost_lambda.*R - Expel(:,3);
uCI(:,4) = uCI_1daypost_top10_lambda.*R - Expel(:,4);

lCI(:,1) = Expel(:,1) - lCI_1daypre_lambda.*R;
lCI(:,2) = Expel(:,2) - lCI_1daypre_top10_lambda.*R;
lCI(:,3) = Expel(:,3) - lCI_1daypost_lambda.*R;
lCI(:,4) = Expel(:,4) - lCI_1daypost_top10_lambda.*R;

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