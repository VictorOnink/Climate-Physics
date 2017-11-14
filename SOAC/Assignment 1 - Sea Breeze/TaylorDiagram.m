% Create a Taylor diagram for the sea breeze model
close all
clear all

lam=linspace(1,16,16);
cor=[0.80836 0.79531 0.80323 0.81178 0.81751 0.82031 0.82054 0.81857 0.81477 0.80944 0.80286 0.79529 0.78693 0.77797 0.76858 0.75889];
CPRMS=[8.34 5.44 3.75 2.78 2.27 2.03 1.95 1.96 2.00 2.06 2.12 2.19 2.25 2.3 2.36 2.4];
STD=[10.8 7.7 5.8 4.6 3.8 3.3 2.8 2.5 2.2 2 1.9 1.7 1.6 1.5 1.4 1.3];

yyaxis left
title('Correlation and CP-RMS of the Sea Breeze Model','FontSize',16)
xlabel('\lambda (\times10^{-4})','FontSize',14)
xlim([0,17])
plot(lam,cor,'*')
ylabel('Correlation','FontSize',14)
yyaxis right
plot(lam,CPRMS,'*')
ylabel('CP-RMS','FontSize',14)

