%This script determines the rate of expelling SARS-CoV-2 via different 
%droplets and aerosols during respiratory activities. It uses Poisson 
%statistics to determine likelihood at different sizes and then multiplies
%with particle profiles from respiratory activities, which were extracted
%from the literature. The delta method was used to determine variance. 

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

v_SARS2_1daypost = 7.00; %SARS-CoV-2 mean viral load, -1 day {log10 copies/ml}
uCI_1daypost = 7.58;
lCI_1daypost = 6.43;
v_SARS2_1daypost_top10 = 11.40318; %SARS-CoV-2 viral load (98th cp), -1 day {log10 copies/ml}
uCI_1daypost_top10 = 12.48296;
lCI_1daypost_top10 = 10.41681;

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
1.093135587
1.17984024
1.273422079
1.374426585
1.483442504
1.601105279
1.728100758
1.865169186
2.013109523
2.1727841
2.345123647
2.531132715
2.731895535
2.948582337
3.182456169
3.434880261
3.70732597
4.001381361
4.275030589
4.567394329
4.954817061
5.293670507
5.742698583
6.135434215
6.555028521
7.111049498
7.636123227
8.283846527
8.940899546
9.455633002
10.20562906
10.95920386
11.82845922
12.63739204
13.63975689
14.572562
15.72841946
16.88979232
18.22944643
19.57549295
21.02093037
22.92033677
24.2398756
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
3.3689E+00
4.6562E+00
6.1692E+00
7.6184E+00
8.8930E+00
9.5403E+00
9.9506E+00
9.4048E+00
8.6422E+00
7.8304E+00
6.8016E+00
5.8253E+00
5.2044E+00
4.3951E+00
4.0960E+00
3.7112E+00
3.2693E+00
2.8001E+00
2.4667E+00
2.1426E+00
1.7346E+00
1.4444E+00
1.1369E+00
9.0759E-01
7.2451E-01
5.3905E-01
4.0107E-01
2.8207E-01
1.9560E-01
1.4760E-01
9.8125E-02
6.4322E-02
4.0993E-02
2.7252E-02
1.7125E-02
1.1069E-02
6.4832E-03
3.9610E-03
2.2556E-03
1.3976E-03
9.4229E-04
6.7209E-04
7.1095E-04
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