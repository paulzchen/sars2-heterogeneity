%This script determines the expected rate of shedding SARS-CoV-2 via 
%droplets and aerosols while breathing. 

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
1.078609041
1.168587651
1.254851023
1.347482227
1.453406513
1.553763084
1.661049197
1.775743331
1.906826017
2.065895288
2.238234271
2.392782578
2.558002324
2.746830186
2.962755873
3.181462029
3.446862668
3.684865811
3.956877031
4.248967706
4.562620074
4.85600163
5.077040429
5.308140622
5.59938819
5.906615958
6.230700729
6.51431402
6.810837014
7.120857344
7.511565228
7.923710508
8.358469416
8.777922362
9.218424709
9.638035261
10.07674593
10.53542613
11.01498486
11.51637246
12.04058257
12.81480545
13.51792942
14.38714735
15.4491848
16.88764783
17.97354241
19.12926122
20.54135435
21.47636974
22.65473725
24.11146211
25.43441213
26.94964454]*1e-6; %hydrated particle diameters for given respiratory activity {m}
V = (p*pi*Vml_Vg/6)* d.^3; %particle volumes {ml}
R = [0.000189985
0.000352872
0.000664703
0.00120033
0.002198296
0.003597262
0.005969944
0.010334891
0.015762654
0.024727358
0.040463504
0.059162907
0.103354972
0.165054716
0.276537901
0.468144459
0.74480181
1.120086805
1.443235714
1.865276003
2.281793084
2.642023869
2.842935297
3.115693828
3.003583118
2.842935297
2.546957027
2.120537811
1.733456166
1.443235714
1.158367999
0.94692023
0.788383993
0.6934815
0.588053505
0.507874451
0.422844581
0.352050667
0.272395128
0.210762577
0.163075104
0.135772579
0.121637283
0.117260459
0.126177474
0.14084038
0.16916199
0.195867898
0.226789916
0.253144851
0.272395128
0.282562455
0.262593648
0.235254998
0.206935948
0.175476072
0.151550518
0.130887131
0.110988738
0.09411544
0.082786202
0.074167322
0.066445755
0.059528081
0.054316789
0.047778352
0.04126394
0.035637745
0.030219841
0.025160344
0.020194163
0.016208213]; %rate of particle being expelled, corresponding to particle diameters {particles/interval}


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

Expel(:,1) = P_SARS2_1daypre*n.*R;
Expel(:,2) = P_SARS2_1daypre_top10* n.*R;
Expel(:,3) = P_SARS2_1daypost* n.*R;
Expel(:,4) = P_SARS2_1daypost_top10* n.*R;

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