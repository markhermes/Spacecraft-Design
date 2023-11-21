
figure(2);
set(gcf,'units','inches','position',[0 0 10 4],'color','w');
subplot(1,2,1);
plot(t,tau(:,1),'.r'); grid on; box on; hold on;
plot(t,tau(:,2),'.g');
plot(t,tau(:,3),'.b');
xlabel('time (s)'); ylabel('Applied Torques (N m)');
title('Torques applied to spacecraft from reaction wheels');
legend({'x-axis','y-axis','z-axis'}); ylim([-0.015 0.015]);
xlim([0 200]);

subplot(1,2,2);
semilogy(t,err); grid on; box on;
xlabel('time (s)'); ylabel('Angular Error (deg)');
title('Angular error between Front-Starboard-Nadir and attitude');
xlim([0 5000]);


