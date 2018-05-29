%% config
config;
meshEdgeColor=[.2,.2,.2];
startColor=[0, 0, 1];
endColor=[1, 0, 0];
fontSize=15;

%% U region
drawURegion=false;
if drawURegion
    %% read data
    load(['resultData/UShaped_ExplictDisabled_' num2str(scale) '.mat']);
    disEaResult=EaResult;
    disEmResult=EmResult;
    disResult=Result;
    load(['resultData/UShaped_' num2str(scale) '.mat']);
    
    %% draw
    hold on;
    patch('faces',t,'vertices',p,'facecolor','w','edgecolor',meshEdgeColor);
    initialLineX=[0.0, 0.4, 0.4, 0.6, 0.6, 1.0];
    initialLineY=[1.0, 1.0, 0.6, 0.6, 1.0, 1.0];
    plot(initialLineX, initialLineY, 'r-', 'LineWidth', 3);
    ax = gca;
    ax.FontSize = fontSize;
    
    pause
    clf;
    hold on;
    xlim([0,1]);
    ylim([0,1]);
    vertices=[0,0; 1,0; 1,1; 0.6,1; 0.6,0.6; 0.4,0.6; 0.4,1; 0,1];
    plot(vertices(:,1), vertices(:,2), 'Color', meshEdgeColor);
    plot([vertices(end,1), vertices(end,1)], [vertices(end,2) vertices(1,2)], 'Color', meshEdgeColor);
    for idxStep=1:12:length(EaResult)
        [~, ha]=tricontour(p, t, AnalyticalResult{idxStep}, [0, 0]);
        for idxH=1:numel(ha)
            ha(idxH).EdgeColor='r';
        end
        
        [~, hd]=tricontour(p, t, disResult{idxStep}, [0, 0]);
        for idxH=1:numel(hd)
            hd(idxH).EdgeColor='b';
        end
        
        [~, hr]=tricontour(p, t, Result{idxStep}, [0, 0]);
        for idxH=1:numel(hr)
            hr(idxH).EdgeColor='k';
        end
    end
    legend([ha(1) hd(1) hr(1)],...
        {'Analytical result',...
        'Simulation result without EC','Simulation result with EC'},...
        'Location','SouthWest', 'FontSize', fontSize);
    ax = gca;
    ax.FontSize = fontSize;
    
    pause;
    clf;
    hold on;
    plot(disEmResult, 'b-.', 'LineWidth', 1);
    plot(disEaResult, 'r-.', 'LineWidth', 1);
    plot(EmResult, 'b-', 'LineWidth', 2);
    plot(EaResult, 'r-', 'LineWidth', 2);
    h=legend('E_m without EC','E_a without EC',...
        'E_m with EC','E_a with EC',...
        'Location','NorthWest');
    h.FontSize = fontSize;
    xlabel('Step','FontSize', fontSize);
    ylabel('Error','FontSize', fontSize);
    ax = gca;
    ax.FontSize = fontSize;

end

%% 3 circles
drawCircles=false;
drawCirclesGap=3;
if drawCircles
    load(['resultData/3Circles_' num2str(scale) '.mat']);
    
    %% iso lines
    idxArray=1:drawCirclesGap:size(Result,1);
    nLines=numel(idxArray);
    
    hold on;
    for idx=1:nLines
        phi=Result{idxArray(idx)};
        [C, h]=tricontour(p, t, phi, [0, 0]);
        crtColor=(startColor * (nLines-idx) + endColor * idx)/nLines;
        for idxH=1:numel(h)
            h(idxH).EdgeColor=crtColor;
        end
    end
    ax = gca;
    ax.FontSize = fontSize;

    
    pause;
    %% error
    clf;
    hold on;
    plot(EmResult);
    plot(EaResult);
    ax = gca;
    ax.FontSize = fontSize;
    h=legend('E_m','E_a');
    h.FontSize = fontSize;
    xlabel('Step','FontSize', fontSize);
    ylabel('Error','FontSize', fontSize);
    ax = gca;
    ax.FontSize = fontSize;
end

%% Rate stick
drawRateStick=false;
drawRateStickGap=8;
if drawRateStick
    load(['resultData/RateStick_' num2str(scaleStick) '.mat']);
    
    %% mesh
    hold on;
    pbaspect([1 0.05 1]);
    patch('faces',t,'vertices',p,'facecolor','w','edgecolor',meshEdgeColor);
    plot(0, 0.025, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    text(-0.05, 0.025, '$\vec{x_0}$', 'Interpreter', 'latex', 'FontSize', fontSize+1);
    ax = gca;
    ax.YAxisLocation = 'right';
    ax.FontSize = fontSize+1;
    
    %% contour
    pause;
    idxArray=1:drawRateStickGap:size(Result,1);
    nLines=numel(idxArray);
    
    clf;
    hold on;
    plot(0, 0.025, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    text(-0.05, 0.025, '$\vec{x_0}$', 'Interpreter', 'latex', 'FontSize', fontSize+1);
    pbaspect([1 0.05 1]);
    for idx=1:nLines
        phi=Result{idxArray(idx)};
        [C, h]=tricontour(p, t, phi, [0, 0]);
        crtColor=(startColor * (nLines-idx) + endColor * idx)/nLines;
        for idxH=1:numel(h)
            h(idxH).EdgeColor=crtColor;
        end
    end
    ax = gca;
    ax.YAxisLocation = 'right';
    ax.FontSize = fontSize+1;

    
    pause;
    %% error
    clf;
    hold on;
    pbaspect auto;
    plot(EmResult);
    plot(EaResult);
    h=legend('E_m','E_a','Location','SouthEast');
    h.FontSize = fontSize;
    xlabel('Step','FontSize', fontSize);
    ylabel('Error','FontSize', fontSize);
    ax = gca;
    ax.FontSize = fontSize;

end

%% Redistance
drawRedistance=true;
if drawRedistance
    x = -1:scale:1;
    y = x;
    [gX, gY] = meshgrid(x, y);
    %gX = reshape(gX, numel(gX), 1);
    %gY = reshape(gY, numel(gY), 1);
    % only draw the z = 0 2D section
    gZ = 0;
    h = 0.5;
    phi = -min(min(min(min(h+gX, h-gX)), min(h+gY, h-gY)), min(h+gZ, h-gZ));
    
    hold on;
    contour(gX, gY, phi, [0, 0], 'r-', 'LineWidth', 2);
    
    phi = phi .* ((gX-0.3).^2 + (gY-0.2).^2 + (gZ-0.1).^2 + 0.005);
    contour(gX, gY, phi, min(min(phi)):0.03:max(max(phi)));
    
    ax = gca;
    ax.FontSize = fontSize;
end
