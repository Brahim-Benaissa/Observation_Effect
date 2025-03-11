clc
clear

% Load data from .mat files
load('Coef.mat');  
load('galaxy_mass_properties');

% Create figure with publication-quality dimensions and white background
fig = figure('Position', [30, 70, 1000, 500], 'Color', 'w');
hold on;

% Define a refined, publication-quality color scheme
bar_face_color = [0.6350, 0.0780, 0.1840];  % Deep, rich red
bar_edge_color = [0.3, 0.0, 0.1];           % Slightly darker red for edges

% Plot the bar chart with refined styling
% 'BarWidth' is adjusted to prevent overlap when many bars are present.
hBar = bar(corr_M_R_curvatures, 'FaceColor', bar_face_color, ...
    'EdgeColor', bar_edge_color, 'LineWidth', 1.5, 'BarWidth', 0.8);

% Enhance axes labels with bold fonts for clarity
xlabel('Galaxy', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Curvature Correlation', 'FontSize', 18, 'FontWeight', 'bold');

% Customize axes for a clean, professional appearance
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'TickDir', 'out', 'Box', 'off');
grid on;
set(gca, 'GridLineStyle', '-', 'GridColor', [0.8, 0.8, 0.8]);

% Adjust x-axis limits for visual balance
xlim([0, numel(corr_M_R_curvatures) + 1]);
ylim([0.7 1]);
% Save the figure with high-quality output 
print(fig, 'Galaxy_Correlation_BarChart.png', '-dpng', '-r600');
% print(fig, 'Galaxy_Correlation_BarChart.jpg', '-djpeg', '-r600');
% savefig(fig, 'Galaxy_Correlation_BarChart.fig');

hold off;
close(fig);



% Create figure with publication-quality dimensions and white background
fig = figure('Position', [30, 70, 1000, 500], 'Color', 'w');
hold on;

% Define a refined, publication-quality color scheme
bar_face_color = [0.6350, 0.0780, 0.1840];  % Deep, rich red
bar_edge_color = [0.3, 0.0, 0.1];           % Slightly darker red for edges

% Plot the bar chart with refined styling
% 'BarWidth' is adjusted to prevent overlap when many bars are present.
hBar = bar(Mezzi_Mass, 'FaceColor', bar_face_color, ...
    'EdgeColor', bar_edge_color, 'LineWidth', 1.5, 'BarWidth', 0.8);

% Enhance axes labels with bold fonts for clarity
xlabel('Galaxy', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Mezzi Mass coefficient', 'FontSize', 18, 'FontWeight', 'bold');

% Customize axes for a clean, professional appearance
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'TickDir', 'out', 'Box', 'off');
grid on;
set(gca, 'GridLineStyle', '-', 'GridColor', [0.8, 0.8, 0.8]);

% Adjust x-axis limits for visual balance
xlim([0, numel(corr_M_R_curvatures) + 1]);
ylim([1 100]);
% Save the figure with high-quality output  
print(fig, 'Mezzi_Mass_coefficient_BarChart.png', '-dpng', '-r600');
% print(fig, 'Mezzi_Mass_coefficient_BarChart.jpg', '-djpeg', '-r600');
% savefig(fig, 'Mezzi_Mass_coefficient_BarChart.fig');

hold off;
close(fig);




for galaxy = 1:length(Mezzi_Mass)

% Create figure with same size and white background
fig = figure('Position', [100 100 800 600], 'Color', 'w');
hold on;

% Replace these with your actual data
X = Mezzi_Scale_curvature{galaxy};       % Your Mezzi/r² data vector
Y = Ricci_curvature{galaxy}; % Your Ricci curvature data vector

% Create scatter plot with styling
scatter(X, Y, 60, 'MarkerEdgeColor', [0.6350, 0.0780, 0.1840],...
    'MarkerFaceColor', [0.9, 0.9, 0.9],...
    'LineWidth', 1.5);

% Calculate linear fit
p = polyfit(X, Y, 1);
y_fit = polyval(p, X);

% Calculate confidence intervals
[~, S] = polyfit(X, Y, 1);
ci = predint(p, X, 0.98);

% Plot fit line and confidence interval
plot(X, y_fit, 'Color', [0.1, 0.2, 0.6], 'LineWidth', 2.5);
fill([X; flipud(X)], [ci(:,1); flipud(ci(:,2))], [0.1, 0.2, 0.6],...
    'FaceAlpha', 0.15, 'EdgeColor', 'none');

% Set axes labels
xlabel('Mezzi Scale Curvature (1/m²)', 'FontSize', 18);
ylabel('Ricci Curvature (1/m²)', 'FontSize', 18);

title(sprintf('Galaxy %d : %s', galaxy, galaxy_names{galaxy}), 'FontSize', 16);

% Set grid and axes properties
grid on;
set(gca, 'GridColor', [0.8, 0.8, 0.8],...
    'FontSize', 16,...
    'LineWidth', 1.5,...
    'Box', 'on');

% Add legend
% legend('Data Points', 'Linear Fit', '98% Confidence Interval','Location', 'north');

% Save the figure
% print(fig, sprintf('Galaxy_%d_Mezzi_R2_vs_Ricci.jpg', galaxy), '-djpeg', '-r600');
print(fig, sprintf('Galaxy_%d_Mezzi_R2_vs_Ricci.png', galaxy), '-dpng', '-r600');


hold off;

close all

end




function ci = predint(p, X, level)
    % Calculate the number of data points
    n = numel(X);
    % Degrees of freedom for linear fit (n - 2)
    df = n - 2;
    
    % Compute predicted values using the linear model
    y_fit = polyval(p, X);
    
    % Retrieve Y from the base workspace (not recommended, but used here for compatibility)
    Y = evalin('base', 'Y');
    
    % Calculate residuals
    residuals = Y - y_fit;
    
    % Estimate the residual standard error (sigma)
    sigma = sqrt(sum(residuals.^2) / df);
    
    % Compute mean of X and sum of squared deviations
    x_bar = mean(X);
    Sxx = sum((X - x_bar).^2);
    
    % Calculate standard error for each prediction
    SE = sigma * sqrt(1/n + (X - x_bar).^2 / Sxx);
    
    % Determine the t-critical value for the given confidence level
    alpha = level;
    t_critical = tinv(1 - alpha/2, df);
    
    % Compute confidence intervals
    delta = t_critical * SE;
    ci = [y_fit - delta, y_fit + delta];
end
