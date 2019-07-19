% plot gravity comparison
%load ../../radialmodels/PREM/prem3L_noocean_gravity.mat

if 0
figure;
plot(rnrm,gnrm,'.');
hold on;
plot(RI,G*gref,'ko');
legend('sol. from FMM','semi-analytic sol.','location','southeast')
xlabel('Radius (km)');
ylabel('Field (m/s^2)');
axis square
end

if 0
% a simple correction
dif = min(gpot)-min(G*phi);
figure
plot(rnrm,gpot-dif,'.');
hold on;
plot(RI,G*phi,'ko');
legend('sol. from FMM','semi-analytic sol.','location','southeast')
xlabel('Radius (km)');
ylabel('Potential (10^3 m^2/s^2)');
axis square
end