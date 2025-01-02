function plot_raster(glmTbl, pressT, freq, sub)


% Extract data for condition 1, condition 2, and condition 5
data_condition1 = glmTbl.(freq)(glmTbl.(freq).condition == 1, :);
data_condition2 = glmTbl.(freq)(glmTbl.(freq).condition == 2, :);
data_condition3 = glmTbl.(freq)(glmTbl.(freq).condition == 5, :);

% Filter rows where the subject column is not equal to the variable 'sub'
data_condition1 = data_condition1(data_condition1.subject == sub, :);
data_condition2 = data_condition2(data_condition2.subject == sub, :);
data_condition3 = data_condition3(data_condition3.subject == sub, :);

% Sort data based on timingPrec in descending order
sorted_condition1 = sortrows(data_condition1, 'timingPrec', 'descend');
sorted_condition2 = sortrows(data_condition2, 'timingPrec', 'descend');
sorted_condition3 = sortrows(data_condition3, 'timingPrec', 'descend');

% Initialize figure
figHandle = figure;

% Define colors using RGB triplets
colors = {[196, 27, 177] / 255, [19, 141, 13] / 255, [249, 100, 14] / 255,[3, 191, 189] / 255}; % Magenta, Green, Orange (RGB), Cyan

% Plot condition 1 data
subplot(3, 1, 1);
hold on;
for i = 1:height(sorted_condition1)
    trial_num = sorted_condition1.trial(i);
    % Retrieve data for this trial from pressT
    press_data = pressT(trial_num, :);
    for j = 1:4 % Iterate over each column
        scatter(press_data(j), i, 36, 'MarkerFaceColor', colors{j}, 'MarkerEdgeColor', colors{j}, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    end
end
xline(800:800:max(pressT(:)), 'k', 'LineWidth', 2); % Add vertical lines every 800ms
% xlabel('Time (s)');
% ylabel('Trials');
% title('Condition 1');
xticks([0, 1000, 2000, 3000]);
xticklabels({''});
xlim([0, 3000]);
ylim([0, 150]);

% Plot condition 2 data
subplot(3, 1, 2);
hold on;
for i = 1:height(sorted_condition2)
    trial_num = sorted_condition2.trial(i);
    % Retrieve data for this trial from pressT
    press_data = pressT(trial_num, :);
    for j = 1:4 % Iterate over each column
        scatter(press_data(j), i, 36, 'MarkerFaceColor', colors{j}, 'MarkerEdgeColor', colors{j}, 'MarkerEdgeColor', 'k', 'LineWidth', 0.2);
    end
end
xline(400:400:1200, 'k', 'LineWidth', 2); % Add vertical lines every 400ms
% xlabel('Time (s)');
% ylabel('Trials');
% title('Condition 2');
xticks([0, 1000, 2000, 3000]);
xticklabels({''});
xlim([0, 3000]);
ylim([0, 150]);

% Plot condition 3 data
subplot(3, 1, 3);
hold on;
for i = 1:height(sorted_condition3)
    trial_num = sorted_condition3.trial(i);
    % Retrieve data for this trial from pressT
    press_data = pressT(trial_num, :);
    for j = 1:4 % Iterate over each column
        scatter(press_data(j), i, 36, 'MarkerFaceColor', colors{j}, 'MarkerEdgeColor', colors{j}, 'MarkerEdgeColor', 'k', 'LineWidth', 0.2);
    end
end
% xlabel('Time (s)');
% ylabel('Trials');
% title('Condition 3');
xticks([0, 1000, 2000, 3000]);
xticklabels({'0', '1', '2', '3'});
xlim([0, 3000]);
ylim([0, 150]);

% Export the figure
exportgraphics(figHandle, sprintf('./figures/raster_sub%i.eps', sub), 'ContentType', 'vector', 'BackgroundColor', 'none');

hold off;

end