%% Format all open figures with plot box aspect ratio

labelFontSize  = 12;
tickFontSize   = 11;
legendFontSize = 10;
lineWidthValue = 1;

% targetColor = [0.47 0.67 0.19];

figs = findall(groot, 'Type', 'figure');

for i = 1:length(figs)

    fig = figs(i);

    if isprop(fig, 'Number')
        figNum = fig.Number;
    else
        figNum = NaN;
    end

    axs = findall(fig, 'Type', 'axes');

    for j = 1:length(axs)

        ax = axs(j);

        if strcmpi(get(ax, 'Tag'), 'legend')
            continue
        end

        % Axis aspect ratio: x:y:z = 5:2:1
        pbaspect(ax, [5 2 1]);

        % Tick label font size
        ax.FontSize = tickFontSize;

        % Axis label font size
        ax.XLabel.FontSize = labelFontSize;
        ax.YLabel.FontSize = labelFontSize;

        if ~isempty(ax.ZLabel)
            ax.ZLabel.FontSize = labelFontSize;
        end

        ax.Title.FontSize = labelFontSize;
        ax.LineWidth = 1;
        ax.Box = 'on';

        % Change line color up to figure 89
        if figNum <= 89
            lines = findall(ax, 'Type', 'line');

            for k = 1:length(lines)
                % lines(k).Color = targetColor;
                lines(k).LineWidth = lineWidthValue;
            end
        end
    end

    % Legend font size
    lgds = findall(fig, 'Type', 'legend');

    for j = 1:length(lgds)
        lgds(j).FontSize = legendFontSize;
    end
end