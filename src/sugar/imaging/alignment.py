# (C) 2024, Tom Eulenfeld, MIT license
from warnings import warn
import matplotlib.pyplot as plt

import numpy as np

from sugar.data import gcode
from sugar.imaging.colors import get_color_scheme
from matplotlib.colors import to_rgb
from matplotlib.patches import Rectangle

def _get_colordict(color, alphabet, default='white', gap='-', gap_color='white'):
    if color is None:
        color = {l: 'C{}'.format(i % 10) for i, l in enumerate(alphabet)}
    try:
        color = to_rgb(color)
    except ValueError:
        if isinstance(color, dict):
            colord = color
        elif isinstance(color, (tuple, list)):
            colord = {l: color[i % len(color)]  for i, l in enumerate(alphabet)}
        else:
            colord = get_color_scheme(color)
    else:
        colord = {l: color for l in alphabet}
    for g in gap:
        colord.setdefault(g, gap_color)
    return {l: to_rgb(colord.get(l, default)) for l in list(alphabet) + list(gap)}


def _get_fts_colordict(fts, fts_color, fts_colorby):
    if isinstance(fts_colorby, str):
        metaget = fts_colorby
        fts_colorby = lambda ft: ft.meta.get(metaget)
    if fts_color is None:
        fts_color = ['C{}'.format(i) for i in range(10)]
    if isinstance(fts_color, str):
        fts_color = [fts_color]
    if not isinstance(fts_color, dict):
        d = {}
        i = 0
        for ft in fts:
            k = fts_colorby(ft)
            if k not in d:
                d[k] = fts_color[i]
                i = (i + 1) % len(fts_color)
        fts_color = d
    for k in fts_color:
        fts_color[k] = to_rgb(fts_color[k])
    return fts_color, fts_colorby


def _get_xy(extent, n, m):
    if extent is None:
        extent = [-0.5, n - 0.5, -0.5, m - 0.5]
    if len(extent) == 4:
        extent = [extent[0:2], extent[2:4]]
    x, y = extent
    if len(x) == 2:
        x = np.linspace(x[0], x[1], n+1)
    if len(y) == 2:
        y = np.linspace(y[1], y[0], m+1)
    return x, y


def _despine(ax, show_spines, spine_offset):
    if show_spines is not True:
        sides = ['top', 'right', 'left', 'bottom']
        if show_spines is False:
            despine = [True, True, True, True]
        elif isinstance(show_spines, str):
            despine = [side == show_spines for side in sides]
        else:
            despine = show_spines
        for despineit, side in zip(despine, sides):
            if despineit:
                ax.spines[side].set_visible(False)
            elif spine_offset:
                ax.spines[side].set_position(('outward', spine_offset))

_FTS_ACCOUNT_FOR_GAPS = True


def plot_alignment(
        seqs, fname=None, *,
        ax=None, figsize=None,
        extent=None, aspect=None,
        gap='- ',
        color='gray', gap_color='white',
        symbols=False,
        symbol_color='black',
        symbol_gap_color='black',
        symbol_size=None,
        scale_symbol_size=1,
        symbol_kw=None,
        fts=None,
        fts_display='facecolor',
        fts_colorby='type',
        fts_color=None,
        fts_color_gap_alpha=1,
        fts_alpha=None,
        fts_box_groups='type',
        fts_box_lw=5,
        fts_box_kw=None,
        show_spines=False, spine_offset=None,
        xticks=True,
        dpi=None, transparent=None, bbox_inches=None, show=False,
        **kw
        ):
    r"""
    Plot an alignment

    :param seqs: sequences
    :param fname: The filename if saving the plot, the figure will be closed afterwards, default: do not save the plot and return figure
    :param ax: The ax to plot in, default: create a new ax
    :param figsize: The size of the created figure, only applies for ``ax=None``
    :param extent: The extent of the plotted alignment in data coordinates ``[xmin, xmax, ymin, ymax]``,
        the default plots each symbol and sequence at integer coordinates
    :param aspect: Wether to shrink the axis to guarantee the given aspect ratio,
        default no shrinkage, ``aspect=2`` gives visually appealing plots if symbols are also plotted.
    :param gap: The characters recognized as gaps, default is ``'- '``
    :param color: The background color,
        might be any constant color (defaults to ``'gray'``),
        a list of colors, or
        None for the default matplotlib color cycle,
        a dictionary mapping the alphabet to colors, or
        a supported color scheme, see `.get_color_scheme()`, e.g. ``'flower'``.
    :param gap_color: The color of gaps, default is ``'white'``
    :param symbols: Wether to plot symbols (letters), default ``False``
    :param symbol_color: The color of the symbols, takes the same values as the color parameter,
        default is ``'black'``
    :param symbol_gap_color: the color of the gap symbol, default is ``'black'``
    :param symbol_size: The font size of the symbols, by default a visually appealing size is calculated
    :param scale_symbol_size: Scale factor for the automatically calculated symbol size, default is 1
    :param symbol_kw: A dictionary of additional parameters passed to matplotlib's :meth:`~matplotlib.axes.Axes.annotate`
    :param fts: Wether to highlight features, defaults to no highlighting,
        might be a FeatureList object or just True to use the features which are attached to the sequences object.
    :param fts_display: How to display the features, one of ``'facecolor'`` (default) and ``'box'``,
        boxes will range over all sequences
    :param fts_colorby: How to define the color of the features, might be any key in the metadata,
        defaults to ``'type'``, but can also be a function taking a Feature and returning an identifier
    :param fts_color: The color of the features, similarly as with the color parameter,
        this might be a constant color,
        a list of colors, or
        None for the default matplotlib color cycle (the default), or
        a dictionary mapping the feature identifiers to colors.
    :param fts_color_gap_alpha:
        The alpha value of the feature color for gaps (default: 1)
    :param fts_alpha: Transparency of the features
    :param fts_box_groups: For the ``fts_display='box'`` option,
       we need to specify which features belong into the same box,
       this parameter can be a list of FeatureList objects,
       alternatively this parameter is passed to `.FeatureList.groupby()` to
       define the groups, defaults to ``'type'``
    :param fts_box_lw: linewidth of the the boxes, default: 5
    :param fts_box_kw: Dictionary of additional parameters passed to matplotlib's `~matplotlib.patches.Rectangle` to create the feature boxes
    :param show_spines,despine_offset: Parameters passed to seaborn's despine function,
        the default ``show_spines=False`` removes axes spines
    :param xticks: True leaves the xticks (default), False turns them off, can also be a list of xticks
    :param dpi,transparent,bbox_inches: Parameters passed to :meth:`~matplotlib.figure.Figure.savefig` if the figure is saved
    :param show: True shows the figure, default: False
    :param \*\*kw: Other kwargs are passed to matplotlib's :meth:`~matplotlib.axes.Axes.pcolormesh`

    :return: Axes object if ``fname=None``, otherwise the figure is saved and closed

    .. rubric:: Example

    >>> from sugar import read
    >>> seqs = read('https://osf.io/download/j2wyv')
    >>> seqs.plot_alignment(show=True, figsize=(10, 4))

    .. image:: ../_static/ali1.png
       :width: 60%

    >>> seqs[:, 70:120].plot_alignment(show=True, color=None, figsize=(10,8),
    ...                                symbols=True, aspect=2, alpha=0.5)

    .. image:: ../_static/ali2.png
       :width: 40%

    >>> seqs2 = seqs[:5, :150].copy()
    >>> seqs2.translate(complete=True).plot_alignment(
    ...     show=True, color='flower', figsize=(10,8),  symbols=True,
    ...     aspect=2, alpha=0.5, edgecolors='w')

    .. image:: ../_static/ali3.png
       :width: 40%
    """
    if gap is None:
        gap = ''
    alphabet = sorted(set(''.join(str(seq) for seq in seqs)) - set(gap))
    color = _get_colordict(color, alphabet, gap=gap, gap_color=gap_color)
    if fts:
        if fts is True:
            fts = seqs.fts
        fts_color, fts_colorby = _get_fts_colordict(fts, fts_color, fts_colorby)
    if not 0 <= fts_color_gap_alpha <= 1:
        raise ValueError('fts_color_gap_alpha has to be a number 0 <= alpha <= 1')
    if fts_color_gap_alpha < 1 and fts_alpha is None:
        fts_alpha = 1
    lens = [len(seq) for seq in seqs]
    n = max(lens)
    if len(set(lens)) > 1:
        warn('fill up short sequences with empty space')
        seqs = seqs.copy().str.ljust(n)

    data = [[color[l] for l in seq.data] for seq in seqs]
    if fts and fts_display == 'facecolor':
        ftsd = fts.groupby('seqid')
        if fts_alpha is not None:
            data_fts = [[(1, 1, 1, 0) for l in seq.data] for seq in seqs]
        for i, seq in enumerate(seqs):
            for ft in ftsd.get(seq.id, []):
                slice = seq.slindex(gap=gap if _FTS_ACCOUNT_FOR_GAPS else None)[ft]
                start, stop, _ = slice.indices(len(data[i]))
                if fts_alpha is None:
                    data[i][slice] = [fts_color[fts_colorby(ft)]] * (stop - start)
                else:
                    data_fts[i][slice] = [fts_color[fts_colorby(ft)] + (fts_alpha,)] * (stop - start)
                if fts_color_gap_alpha < 1:
                    assert fts_alpha is not None
                    for si in range(start, stop):
                        if seq.data[si] in gap:
                            data_fts[i][si] = fts_color[fts_colorby(ft)] + (fts_color_gap_alpha,)
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0, 0, 1, 1])
    else:
        fig = ax.get_figure()
    x, y = _get_xy(extent, n, len(data))
    ax.pcolormesh(x, y, np.array(data), **kw)
    if fts and fts_alpha is not None and fts_display == 'facecolor':
        ax.pcolormesh(x, y, np.array(data_fts), **kw)
    if fts and fts_display == 'box':
        if fts_box_kw is None:
            fts_box_kw = {}
        fts_box_kw.setdefault('fill', False)
        if not isinstance(fts_box_groups, (tuple, list)):
            fts_box_groups = fts.groupby(fts_box_groups).values()
        for ftgroup in fts_box_groups:
            bx1, bx2 = ftgroup.loc_range
            ax.add_patch(Rectangle((x[bx1], min(y)), x[bx2] - x[bx1], abs(y[-1] - y[0]),
                                   color=fts_color[fts_colorby(ftgroup[0])], alpha=fts_alpha,
                                   lw=fts_box_lw, **fts_box_kw))
    if aspect is not None:
        aspect = abs(aspect * len(data) / n / (y[-1] - y[0]) * (x[-1] - x[0]))
        ax.set_aspect(aspect)
    if symbols:
        if symbol_kw is None:
            symbol_kw = {}
        symbol_color = _get_colordict(symbol_color, alphabet, default='black', gap=gap, gap_color=symbol_gap_color)
        symbol_kw.setdefault('family', 'monospace')
        if 'verticalalignment' not in symbol_kw:
            symbol_kw.setdefault('va', 'center_baseline')
        if 'horizontalalignment' not in symbol_kw:
            symbol_kw.setdefault('ha', 'center')
        if symbol_size is None:
            fig.draw_without_rendering()  # to set automatic axis limits
            # transform bounding box of ax to inch
            bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            ratio_covered = abs((x[-1] - x[0]) / (ax.get_xlim()[1] - ax.get_xlim()[0]))
            # Calculate the width of a single symbol box in inch.
            # Font size is the font height specified in pixels with PPI (pixels per inch) of 72.
            # With a typical aspect ratio of monospace fonts of 5:3, we would need a factor of 120 to fill the width of the symbol box
            # We take slightly lower value of 100 here.
            symbol_size = bbox.width * ratio_covered / n * 100 * scale_symbol_size
        elif scale_symbol_size != 1:
            warn('scale_symbol_size parameter is ignored if symbol_size is specified')
        for i in range(len(data)):
            for j in range(n):
                xy = 0.5 * (x[j] + x[j+1]), 0.5 * (y[i] + y[i+1])
                l = seqs[i].data[j]
                ax.annotate(l, xy, color=symbol_color[l], size=symbol_size, **symbol_kw)
    _despine(ax, show_spines, spine_offset)
    ax.set_yticks([])
    if xticks is not True:
        if xticks is False:
            xticks = []
        ax.set_xticks(xticks)
    if fname is not None:
        fig.savefig(fname, dpi=dpi, transparent=transparent, bbox_inches=bbox_inches)
        if show:
            plt.show()
        plt.close(fig)
    else:
        if show:
            plt.show()
        return ax


if __name__ == '__main__':
    from sugar import read

    seqs = read().sl(update_fts=True)[:, :100]
    seqs[1][:10] = '-' * 10
    seqs[1][:5] = ' ' * 5
    print(seqs)
    seqs.fts = seqs.fts.slice(10, 15)[:1]
    print(seqs.fts)
    plot_alignment(seqs, fts=True, aspect=2, rasterized=True, color='0.8', symbols=True, symbol_color='flower')
