%This **empty** script calculates the the total rate of shedding SARS-CoV-2
%(i.e., cumulative rate across all particle sizes) via a respiratory 
%activity across case percentiles during the infectious period. Data for 
%the percentiles and 95% CIs were determined based on Weibull quantile 
%functions for respiratory viral load on each day and likelihood modelling.

%As an empty script, the raw data can be input into the arrays below to
%solve for the total shedding rates.

clear all
clc
close all

%Define values
gamma = 0.1e-2; %viability proportion {%}
x = [1:1:99]*1e-2; %proportions from 1% to 99%


%{
**empty arrays to input raw data (see uploaded spreadsheets)**
R_vol = ; %Total volumetric rate of particles expelled {m^3/min}
vl = []; %estimate of respiratory viral load {log10 copies/ml}     
              %from 1% to 99%
lCI_vl = []; %lower 95% CI for estimate {log10 copies/ml}
uCI_vl = []; %upper 95% CI for estimate {log10 copies/ml}
%}




%Initialize matrices for respiratory viral load, CIs and Poisson parameters
Expel = zeros(length(vl),1); %Poisson means
lCI_lambda = zeros(length(vl),1); %lower 95% CI Poisson means
uCI_lambda = zeros(length(vl),1); %upper 95% CI Poisson means

%Solve for respiratory viral loads
for i = 1:length(vl)
    Expel(i) = (10^vl(i))* R_vol*gamma;
    lCI_lambda(i) = (10^lCI_vl(i))* R_vol*gamma;
    uCI_lambda(i) = (10^uCI_vl(i))* R_vol*gamma;
end

del_lCI = Expel - lCI_lambda; %lower 95% CI difference for each expel
del_uCI = uCI_lambda - Expel; %upper 95% CI difference for each expel

Output(:,1) = Expel;
Output(:,2) = del_lCI;
Output(:,3) = del_uCI;


figure
plot(x,Expel)
hold on
plot(x,Expel + del_uCI)
plot(x,Expel - del_lCI)
hold off
title('SARS-CoV-2 total shedding rate vs percentile')
ylabel('Total shedding rate (virions/unit)')
xlabel('Percentile (%)')
set(gca,'YScale','log')