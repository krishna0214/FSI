function graph_V(x, u_star)
    figure(1);
    % Plot the graph
    plot(x, u_star, 'LineWidth', 2); % Plot with thicker line
    hold on; % Keep the plot for adding more plots or annotations

    % Add labels and title
    xlabel('x'); % X-axis label
    ylabel('velocity'); % Y-axis label
    title('1. velocity vs x'); % Title
    grid on;

    % Pause for a while
    pause(0.5); % Pause for 5 seconds
    
end


