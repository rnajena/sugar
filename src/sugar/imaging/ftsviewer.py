import sugar
from sugar.core.util import _pop_kws_for_func


def plot_ftsviewer(
        fts: sugar.FeatureList, fname=None, *,
        seqs=None, groupby='seqid',
        figsize=None, ncols=1, sharex=False, sharey=False, fig_kw={}, ax=None,
        axlabel=True, axlabel_kw={},
        same_colors=True, colorby='type', color=None,
        dpi=None, transparent=None, bbox_inches=None, show=False,
        **kw):
    """
    Plot features using DNAFeaturesViewer_ and `~.FeatureList.toftsviewer()`

    :param fts: features
    :param fname: The filename if saving the plot, the figure will be closed afterwards, default: do not save the plot and return figure
    :param seqs: The corresponding sequences, instance of `.BioBasket`,
        will be used to extract the length of each sequence,
        they are also needed when the ``plot_sequence`` option is used.
    :param groupby: Features are grouped by this key, see `~.FeatureList.groupby()`,
        and each group is plotted in a separate ax, defaults to ``'seqid'``.
    :param figsize: The size of the created figure,
        by default a reasonable value is calculated with ``(6 * ncols, 2 * nrows)``,
        has to be tuned for non-usual cases.
    :param ncols: Number of columns in the figure,
        by default all groups are plotted in a single column
    :param sharex,sharey: parameters passed to `plt.subplots <matplotlib.pyplot.subplots>`
    :param fig_kw: dict of other parameters passed to `plt.subplots <matplotlib.pyplot.subplots>`
    :param ax: Instead of specifying the figure with figsize, etc. one may pass an axes instance to draw
        the annotation in, make sure to only have a single annotation group
    :param axlabel: Plot the name of each group in the axis, default: True
    :param axlabel_kw: Dict of corresponding options passed to :meth:`~matplotlib.axes.Axes.annotate`
    :param same_colors: Wether to use the same color mapping in all subplots, default True.
    :param dpi,transparent,bbox_inches: Parameters passed to :meth:`~matplotlib.figure.Figure.savefig` if the figure is saved
    :param show: True shows the figure, default: False
    :param \*\*kw: Other kwargs are passed to
        `~.FeatureList.toftsviewer()`,
        `~dna_features_viewer.GraphicFeature` or
        `~dna_features_viewer.GraphicRecord` or
        `~dna_features_viewer.CircularGraphicRecord`, and
        `GraphicRecord.plot() <dna_features_viewer.GraphicRecord.plot>`.

    :return: Figure object if ``fname=None``, otherwise the figure is saved and closed

    .. rubric:: Example

    >>> from sugar import read
    >>> seqs = read()
    >>> seqs.plot_ftsviewer()

    .. image:: ../_static/fts1.png
       :width: 40%

    """
    try:
        from dna_features_viewer import GraphicRecord
    except ImportError as ex:
        raise ImportError('Please install dna_features_viewer to use ftsviewer functionality') from ex
    import matplotlib.pyplot as plt
    grouped = fts.groupby(groupby, flatten=True)
    N = len(grouped)
    if N > 1 and same_colors:
        from sugar.imaging.alignment import _get_fts_colordict
        color, colorby = _get_fts_colordict(fts, color, colorby)
    if ax is None:
        nrows = (N + ncols - 1) // ncols
        if figsize is None:
            figsize = (6 * ncols, 2 * nrows)
        fig, axes = plt.subplots(nrows, ncols, sharex=sharex, sharey=sharey, squeeze=False, figsize=figsize, **fig_kw)
    else:
        if N > 1:
            raise ValueError('Use ax only with a single annotation group, please')
        fig = ax.figure
        axes = None
    seq = None
    kw2 = _pop_kws_for_func(kw, GraphicRecord.plot)
    for i, gkey in enumerate(grouped):
        if axes is not None:
            ax = axes[i % nrows, i // nrows]
        if seqs is not None:
            seq = seqs.d[grouped[gkey][0].seqid]
        record = grouped[gkey].toftsviewer(seq=seq, colorby=colorby, color=color, **kw)
        record.plot(ax=ax, **kw2)
        if axlabel and gkey is not None:
            axlabel_kw = axlabel_kw.copy()
            axlabel_kw.setdefault('xy', (0, 0.8))
            axlabel_kw.setdefault('xycoords', 'axes fraction')
            axlabel_kw.setdefault('va', 'top')
            ax.annotate(' '.join(map(str, gkey)), **axlabel_kw)
    if axes is not None:
        for i in range(i+1, ncols*nrows):
            axes[i % nrows, i // nrows].axis('off')
    if fname is not None:
        fig.savefig(fname, dpi=dpi, transparent=transparent, bbox_inches=bbox_inches)
        if show:
            plt.show()
        plt.close(fig)
    else:
        if show:
            plt.show()
        return fig


if __name__ == '__main__':
    from sugar import read, read_fts
    fts = read_fts()
    plot_ftsviewer(fts)
    plot_ftsviewer(fts, groupby='type')
    seqs = read()
    plot_ftsviewer(seqs.fts, groupby=('seqid', 'type'))
    seqs.plot_ftsviewer()
