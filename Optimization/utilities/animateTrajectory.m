function animateTrajectory(tspan, x, p, dt)
    keypoints_fn = casadi.Function.load('codegen/keypoints_fn.casadi');
    % Prepare plot handles
    hold on
    h_OB = plot([0],[0],'LineWidth',2);
    h_AC = plot([0],[0],'LineWidth',2);
    h_BD = plot([0],[0],'LineWidth',2);
    h_CE = plot([0],[0],'LineWidth',2);
    
    
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
    
        rA = keypoints(:,1); % Vector to base of cart
        rB = keypoints(:,2);
        rC = keypoints(:,3); % Vector to tip of pendulum
        rD = keypoints(:,4);
        rE = keypoints(:,5);
    
        set(h_title,'String',  sprintf('t=%.2f',t) ); % update title
        
        set(h_OB,'XData',[0 rB(1)]);
        set(h_OB,'YData',[0 rB(2)]);
        
        set(h_AC,'XData',[rA(1) rC(1)]);
        set(h_AC,'YData',[rA(2) rC(2)]);
        
        set(h_BD,'XData',[rB(1) rD(1)]);
        set(h_BD,'YData',[rB(2) rD(2)]);
        
        set(h_CE,'XData',[rC(1) rE(1)]);
        set(h_CE,'YData',[rC(2) rE(2)]);
    
        pause(dt)
    end
end
    