function tight_layout(fsize, font)
% Tight the layout of the figure.
%   fsize (optional): the size of the output figure.
%   font (optional): the size of the output font.
if exist('fsize', 'var')
    set(gcf, 'units', 'normalized', 'outerposition', ...
        [(1-fsize(1))/2, (1-fsize(2))/2, fsize(1) fsize(2)]);
end
get_axes = findall(gcf, 'type', 'axes')';
for ax = get_axes
    if exist('fsize', 'var')
        set(ax, 'fontsize', font);
    end
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    ts = 0.05;
    ti = [ti(1) + ts, ti(2) + ts, ti(3) + ts, ti(4) + ts];
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end
end