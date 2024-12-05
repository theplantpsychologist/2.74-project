
% 
dth = deg2rad(4);
dang = 2;
slopevec = deg2rad(10):dth:deg2rad(50);
angvec = 5:dang:15;
Data = sweeps(slopevec, angvec);
energy_per_height = Data.energy./Data.y;
slopevecdeg = slopevec./pi*180;
close all;
% figure;hold on
heatmap(angvec,flip(slopevecdeg),flipud(energy_per_height))

% xlabel('W (rad/s)'); ylabel('Slope (deg)');
% plot(Slope_val_deg,energy_per_height, "ro-")
% plot(Slope_val_deg,energy, "bo-")
% legend()
% for k = 1:size(energy_per_height,1)
%     hold on
%     plot(angvec, flipud(energy_per_height(k,:)),'o-');
% end
% xlabel('W (rad/s)'); ylabel('Work per unit height (J/m)');
% xlim([5 15])
% ylim([min(min(energy_per_height))-20 max(max(energy_per_height))+20])

function Data = sweeps(slope_vec, ang_vec)
    n = length(slope_vec);
    m = length(ang_vec);
    y = zeros(n,m); energy = zeros(n,m); slope = zeros(n,m); w = zeros(n,m);

    for j = 1:m
        for i = 1:n
            output = simulate_sweep_two_legs(slope_vec(i), ang_vec(j));
            y(i,j) = output(1);
            energy(i,j) = output(2);
            % slope(i,j) = output(3);
            % w(i,j) = output(4);
        end
    end
    Data = table(y, energy);
end

% function Data = sweeps(slope_vec, ang_vec)
%     n = length(slope_vec);
%     m = length(ang_vec);
%     y = zeros(n,m); energy = zeros(n,m); slope = zeros(n,m); w = zeros(n,m);
% 
%     for j = 1:m
%         for i = 1:n
%             output = simulate_sweep_two_legs(slope_vec(i), ang_vec(j));
%             y(i,j) = output(1);
%             energy(i,j) = output(2);
%             % slope(i,j) = output(3);
%             % w(i,j) = output(4);
%         end
%     end
%     Data = table(y, energy);
% end