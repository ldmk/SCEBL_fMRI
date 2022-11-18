function [stats, inds]  = colorpie_stats_lk5_general(v1, v2, name, label1, label2)
% function [stats, inds]  = colorpie_stats_lk(var1, var2, name, label1, label2)
% assesses similarity and differences in bivariate distribution of
% normalized variables var1 (plotted on Y axis) and var2 (plotted on X)
% stats to output: av distance from zero; number of points in each sector,
% indices for each sector (so they can be used for plotting later)...
% name = name or title of figures
% count voxels etc in each quadrant
% first look at voxels that hav weights for CS > 0 (angles between 270-90
% degrees), otherwise tand is not specific
%
% Use:
% [stats, inds] = colorpie_stats_lk5_drug_food(data1.dat, data2.dat, 'Whole brain', 'label1', 'label2');
%
% by Leonie Koban, 2019
% Please refer to the following paper when using this script:
% Koban, L., Jepma, M., López-Solà, M., & Wager, T. D. (2019). Different brain networks mediate the effects of social and conditioned expectations on pain. Nature communications, 10(1), 1-13.

%% figure to illustrate sectors

if size(v1)~=size(v2)
    error('Size of variables do not match')
end

upperh = find(v1 > 0); % find 'upper half' voxels
ta(upperh) = (v2(upperh)./v1(upperh));

inds.sector1u = upperh(find(ta(upperh)>tand(337.5) & ta(upperh)<tand(22.5)));

inds.sector2u = upperh(find(ta(upperh)>tand(22.5) & ta(upperh)<tand(67.5)));

inds.sector3u = upperh(find(ta(upperh)>tand(67.5) & ta(upperh)<tand(90)));

inds.sector7u = upperh(find(ta(upperh)>tand(270) & ta(upperh)<tand(292.5)));

inds.sector8u = upperh(find(ta(upperh)>tand(292.5) & ta(upperh)<tand(337.5)));

lowerh = find(v1 < 0); % find 'lower half' voxels
ta(lowerh) = (v2(lowerh)./v1(lowerh));

inds.sector3l = lowerh(find(ta(lowerh)>-Inf & ta(lowerh)<tand(112.5)));

inds.sector4l = lowerh(find(ta(lowerh)>tand(112.5) & ta(lowerh)<tand(157.5)));

inds.sector5l = lowerh(find(ta(lowerh)>tand(157.5) & ta(lowerh)<tand(202.5)));

inds.sector6l = lowerh(find(ta(lowerh)>tand(202.5) & ta(lowerh)<tand(247.5)));

inds.sector7l = lowerh(find(ta(lowerh)>tand(247.5) & ta(lowerh)< Inf));

inds.sector3c = [inds.sector3u; inds.sector3l];
inds.sector7c = [inds.sector7u; inds.sector7l];


%% summarize indices based on functional relationship

inds.unique_var1 = [inds.sector1u; inds.sector5l];
inds.unique_var2 = [inds.sector3c; inds.sector7c];
inds.shared_pos = [inds.sector2u; inds.sector6l];
inds.shared_neg = [inds.sector8u; inds.sector4l];

%% average distance from origin in each sector

secn = {'sector1u'; 'sector2u'; 'sector3c'; 'sector4l'; 'sector5l'; 'sector6l'; 'sector7c'; 'sector8u'};
stats.names = secn;
%     stats.avdistallsectors = nanmean(sqrt(var1.^2 + var2.^2));
%     d = stats.avdistallsectors;

for sc = 1:8
    stats.numel_total(1,sc) = numel(inds.(secn{sc})); % number of voxels in sector
    stats.proport_total(1,sc) = numel(inds.(secn{sc}))/numel(v1); % proportion of voxels in sector (compared to total number of voxels)
    stats.dist0.(secn{sc}) = sqrt(v1(inds.(secn{sc})).^2 + v2(inds.(secn{sc})).^2); % voxel-wise distance from origin per octant
    stats.dist_av(1,sc) = nanmean(stats.dist0.(secn{sc})); % average distance from origin by octant
%         stats.numel_gtd(1,sc) = numel(find(stats.dist0.(secn{sc}) > d)); % number of voxels that are more than 3 STD from origin
    stats.sumsquareddist0(1,sc) = sum(v1(inds.(secn{sc})).^2 + v2(inds.(secn{sc})).^2);
end

%% colorwheel plot

cols = {[.4 .8 .8];  [.8 .7 .4]; [.8 .3 .1];[.8 .8 .8]; ...
        [.6 1 1]; [1 .9 .5]; [1 .5 .2]; [.8 .8 .8]}; % pink/violet for dim 1, green for dim 2, blue pos diagonal, gray for neg diagonal

scf = mean(stats.dist_av(1,:)) * 10000; % scaling factor

%     scatter(var2(inds.unique_var2), var1(inds.unique_var2), scf*sqrt(var2(inds.unique_var2).^2 + var1(inds.unique_var2).^2), '.', 'MarkerEdgeColor', cols{3}); hold on
%     scatter(var2(inds.unique_var1), var1(inds.unique_var1), scf*sqrt(var2(inds.unique_var1).^2 + var1(inds.unique_var1).^2), '.', 'MarkerEdgeColor', cols{1}); hold on
%     scatter(var2(inds.shared_pos), var1(inds.shared_pos), scf*sqrt(var2(inds.shared_pos).^2 + var1(inds.shared_pos).^2), '.', 'MarkerEdgeColor', cols{2}); hold on
%     scatter(var2(inds.shared_neg), var1(inds.shared_neg), scf*sqrt(var2(inds.shared_neg).^2 + var1(inds.shared_neg).^2), '.', 'MarkerEdgeColor', cols{4}); hold on

figure('Name', ['ColorWheel scatter plots', name])

hp = subplot(1,2,1);
scatter(v2(inds.sector1u), v1(inds.sector1u), 2, '.', 'MarkerEdgeColor', cols{1}); hold on
scatter(v2(inds.sector2u), v1(inds.sector2u), 2, '.', 'MarkerEdgeColor', cols{2}); hold on
scatter(v2(inds.sector3c), v1(inds.sector3c), 2, '.', 'MarkerEdgeColor', cols{3}); hold on
scatter(v2(inds.sector4l), v1(inds.sector4l), 2, '.', 'MarkerEdgeColor', cols{4}); hold on
scatter(v2(inds.sector5l), v1(inds.sector5l), 2, '.', 'MarkerEdgeColor', cols{5}); hold on
scatter(v2(inds.sector6l), v1(inds.sector6l), 2, '.', 'MarkerEdgeColor', cols{6}); hold on
scatter(v2(inds.sector7c), v1(inds.sector7c), 2, '.', 'MarkerEdgeColor', cols{7}); hold on
scatter(v2(inds.sector8u), v1(inds.sector8u), 2, '.', 'MarkerEdgeColor', cols{8}); hold on


sf = .8 .* max(abs([v1;v2])); % to scale circles appropriately

c1 = circle([0 0], .25*sf); hold on
c1.Color = [.6 .6 .6];
c2 = circle([0 0], .5*sf); hold on
c2.Color = [.5 .5 .5];
c3 = circle([0 0], .75*sf); hold on
c3.Color = [.6 .6 .6];
c4 = circle([0 0], sf); hold on
c4.Color = [.5 .5 .5];

xlim([min(v1)*1.2 max(v1)*1.2]); ylim([min(v2)*1.2 max(v2)*1.2]);

xlabel(label2); ylabel(label1);

set(hp, 'FontSize', 14, 'Units', 'points', 'Position', [60, 40, 150 150])



%% figures

ma = .5; % marker alpha
ms = 0.001; % marker size

subplot(1,2,2); 
barplot_columns(stats.sumsquareddist0, 'nofig', 'nostars', 'noviolin', 'MarkerAlpha', ma, 'MarkerSize', ms, 'within', 'color', cols, 'LineWidth', 1.5);
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'Units', 'points', 'Position', [275, 40, 150, 150]); xlabel(''); ylabel('SSD from origin')

% set(gcf, 'Position', [100 100 480 200])


end