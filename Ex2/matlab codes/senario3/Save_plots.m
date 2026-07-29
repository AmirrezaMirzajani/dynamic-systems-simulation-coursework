% Save all open figures as .fig files
figs = findobj('Type', 'figure');
for i = 1:length(figs)
    savefig(figs(i), sprintf('figure_%d.fig', i));
end

% Save all open figures as PNG files with 300 DPI
% figs = findobj('Type', 'figure');
% for i = 1:length(figs)
%     exportgraphics(figs(i), sprintf('figure_%d.png', i), 'Resolution', 300);
% end


figs = findobj('Type', 'figure');
for i = 1:length(figs)
    name = figs(i).Name;
    if isempty(name)
        name = sprintf('figure_%d', i);
    else
        % replace invalid characters for filename
        name = strrep(name, ' ', '_');
        name = strrep(name, '(', '');
        name = strrep(name, ')', '');
        name = strrep(name, '/', '_');
    end
    exportgraphics(figs(i), [name, '.png'], 'Resolution', 300);
end