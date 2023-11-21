for j = 1:3  
    figure(1); subplot(1,3,j);
    plot(time,ant_gyro(:,j)*180/pi,'-g');hold on; box on; grid on;
    plot(time,ant_gyro_removed_bias(:,j)*180/pi,'-k'); 
    plot(time,ant_gyro_true(:,j)*180/pi,'-r');
    if(j==1) ylabel('angular rates (deg/s)'); end
    xlabel('time (s)');
end
subplot(1,3,2);
title('Kalman Filter Gyroscope Bias-Corrected Measurements');
subplot(1,3,3);
legend({'meas gyro rate','true body rate','est body rate'})
set(gcf,'units','inches','position',[0 0 6 3],'color','w');


figure(2);
semilogy(time(2:end), uncorrected_euler_err(2:end)); hold on;
semilogy(time(2:end),kf_attitude_error(2:end)); 
title('Kalman Filter Estimated Attitude Error');
legend({'uncorrected euler error','kf corrected error'})
ylabel('RMS error of euler rates'); xlabel('time (s)');
set(gcf,'units','inches','position',[0 0 6 3],'color','w');
grid on;box on; yticks(flip([1e3 1e2 1e1 1 1e-1 1e-2 1e-3 1e-4 1e-5]))
ylim([1e-5 1e3]);