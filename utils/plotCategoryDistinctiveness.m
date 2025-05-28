function [] = plotCategoryDistinctiveness(ROI, xlabs, data, errors)
hBar = bar(data, 'facecolor', 'flat');
hold on

ylim([-0.03 0.35]);

colors = [0.2 0.647 0.9; 0 0.447 0.741; 1 0.425 0.198; 0.85 0.325 0.098];
for i = 1:8
    for c = 1:length(colors)
        hBar(c).CData(i,:) = colors(c,:);
    end
end

ngroups = size(data, 1);
nbars = size(data, 2);
groupwidth = min(0.8, nbars / (nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i - 1) * groupwidth / (2 * nbars);
    er = errorbar(x, data(:, i), errors(:, i), '.');
    er.Color = [0 0 0];
end

if strcmp(ROI, 'IF_DO')
    ROI = 'VLPFC';
end
title(sprintf('%s: Distinctiveness vs. Category\n', ...
    replace(ROI, '_', '\_')), 'FontSize', 16);
ax = hBar.Parent;
set(ax, 'xtick', 1:length(data));
set(gca, 'xticklabel', xlabs, 'FontSize', 16);
set(gcf, 'position', [100, 100, 1000, 700]);
xlabel(sprintf('\nCategory'), 'FontWeight', 'bold', 'FontSize', 16);
ylabel(sprintf('Distinctiveness\n'), 'FontWeight', 'bold', 'FontSize', 16);
legend({'dynamic', '', 'static', ''}, 'Location', 'northwest');

filename = [ROI '_cat_distinctiveness_subplot.fig'];
saveas(gcf, filename);
hold off
close
end
