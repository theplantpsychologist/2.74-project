function animateTrajectory(tspan, x, p, dt)
    keypoints_fn = casadi.Function.load('codegen/keypoints_fn.casadi');
    % Prepare plot handles
    hold on
    h_OB = plot([0],[0],'LineWidth',2);
    h_AC = plot([0],[0],'LineWidth',2);
    h_BD = plot([0],[0],'LineWidth',2);
    h_CE = plot([0],[0],'LineWidth',2);
    h_OB2 = plot([0],[0],'LineWidth',2);
    h_AC2 = plot([0],[0],'LineWidth',2);
    h_BD2 = plot([0],[0],'LineWidth',2);
    h_CE2 = plot([0],[0],'LineWidth',2);
    
    
    xlabel('x'); ylabel('y');
    h_title = subtitle('t=0.0s');
    
    axis equal
    axis([-.2 .2 -.3 .1]);
    
    % Step through and update animation
    for i = 1:length(tspan)
        % % skip frame.
        % if mod(i, 1)
        %     continue;
        % end
        t = tspan(i);
        z = x(:,i); 
        keypoints = full(keypoints_fn(z,p));
    
        rO = keypoints(:,1); % Vector to base of cart
        rA = keypoints(:,2);
        rB = keypoints(:,3); % Vector to tip of pendulum
        rC = keypoints(:,4);
        rD = keypoints(:,5);
        rE = keypoints(:,6); % Vector to base of cart
        rA2 = keypoints(:,7);
        rB2 = keypoints(:,8); % Vector to tip of pendulum
        rC2 = keypoints(:,9);
        rD2 = keypoints(:,10);
        rE2 = keypoints(:,11);
    
        set(h_title,'String',  sprintf('t=%.2f',t) ); % update title
        
        set(h_OB,'XData',[rO(1) rB(1)]);
        set(h_OB,'YData',[rO(2) rB(2)]);
        
        set(h_AC,'XData',[rA(1) rC(1)]);
        set(h_AC,'YData',[rA(2) rC(2)]);
        
        set(h_BD,'XData',[rB(1) rD(1)]);
        set(h_BD,'YData',[rB(2) rD(2)]);
        
        set(h_CE,'XData',[rC(1) rE(1)]);
        set(h_CE,'YData',[rC(2) rE(2)]);

        set(h_OB2,'XData',[rO(1) rB2(1)]);
        set(h_OB2,'YData',[rO(2) rB2(2)]);
        
        set(h_AC2,'XData',[rA2(1) rC2(1)]);
        set(h_AC2,'YData',[rA2(2) rC2(2)]);
        
        set(h_BD2,'XData',[rB2(1) rD2(1)]);
        set(h_BD2,'YData',[rB2(2) rD2(2)]);
        
        set(h_CE2,'XData',[rC2(1) rE2(1)]);
        set(h_CE2,'YData',[rC2(2) rE2(2)]);
    
        pause(dt)
    end
end
    