% plot gravity comparison
load ../deal_prem/prem3L_noocean_gravity.mat

if 1
figure;
plot(rnrm,gnrm,'.');
hold on;
plot(RI,G*gref,'ko');
legend('sol. from FMM','semi-analytic sol.','location','southeast')
xlabel('Radius (km)');
ylabel('Field (m/s^2)');
axis square
end

if 1
figure
plot(rnrm,gpot,'.');
hold on;
plot(RI,G*phi,'ko');
legend('sol. from FMM','semi-analytic sol.','location','southeast')
xlabel('Radius (km)');
ylabel('Potential (10^3 m^2/s^2)');
axis square
end