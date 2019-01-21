# Set the prefix used for figure labels
fig_label_prefix = 'fig'
# Track figure numbers to create unique auto-generated names
fig_count = 0

def savefig(name='', legend=False, fig=None, ext='.pdf'):
    '''
    Save the current figure (or `fig`) to file using `plt.savefig()`.
    If called with no arguments, automatically generate a unique filename.
    Return the filename.
    '''
    # Get name (without extension) and extension
    if not name:
        global fig_count
        # `pytex.id` is a unique identifier for the session;
        # it's `<family>_<session>_<restart>`
        name = 'auto_fig_{}-{}'.format(pytex.id, fig_count)
        fig_count += 1
    else:
        if len(name) > 4 and name[:-4] in ['.pdf', '.svg', '.png', '.jpg']:
            name, ext = name.rsplit('.', 1)

    # Get current figure if figure isn't specified
    if not fig:
        fig = gcf()

    # Put a nice legend on top of the figure (you may need to adjust this!)
    # Only create legend if axis has labels
    if legend and len(fig.gca().get_legend_handles_labels()[0]) != 0: 
        for ax in fig.axes:
            if ax.is_first_row():
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width, box.height*0.9])
        leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.04), ncol=3, bbox_transform=fig.transFigure, frameon=False)
    fig.savefig(name + ext)
    fig.clf()
    return name

def latex_environment(name, content='', option=''):
    '''
    Simple helper function to write the `\begin...\end` LaTeX block.
    '''
    return '\\\\begin{%s}%s\n%s\n\\\\end{%s}' % (name, option, content, name)

def latex_figure(name=None, caption='', label='', width=0.8):
    ''''
    Auto wrap `name` in a LaTeX figure environment.
    Width is a fraction of `\textwidth`.
    '''
    if not name:
        name = savefig()
    content = '\\\\centering\n'
    content += '\\\\includegraphics[width=%f\\\\textwidth]{%s}\n' % (width, name)
    if not label:
        label = name
    if caption and not caption.rstrip().endswith('.'):
        caption += '.'
    # `\label` needs to be in `\caption` to avoid issues in some cases
    content += "\\\\caption{%s\\\\label{%s:%s}}\n" % (caption, fig_label_prefix, label)
    return latex_environment('figure', content, '[htp]')