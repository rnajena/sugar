from warnings import warn

import sugar
from sugar.core.util import _pop_kws_for_func


def _align_fts(grouped, align, *,
               crop=False,
               inplace=True,
               return_details=False
               ):
    """
    Align features in each group by shifting their locations.

    :param grouped: Grouped features, as returned by `~.FeatureList.groupby()`
    :param align: Align features by this method or list of features,
        can be one of the strings ``'center'``, ``'first'``, ``'last'``,
        or a list or BioBasket of features to align each group to the the first found feature in the list.
    :param crop: Crop the sequence lines to the first and last feature in each group after alignment,
        can be ``True`` to crop to the minimal region including all features,
        or an integer to add some extra space before the first features.
    :param inplace: If True, modify the grouped features inplace, otherwise make a copy first.
    :param return_details: If True, return a tuple of (grouped, seq_start, pos_first_ft, stop_pos_last_ft)
        where seq_start is a dict mapping group keys to the new sequence start (after alignment),
        and pos_first_ft is a dict mapping group keys to the position of the first feature (after alignment),
        and stop_pos_last_ft is a dict mapping group keys to the stop position of the last feature (after alignment).
    :return: Aligned and grouped features, or a tuple if ``return_details=True``
    """
    if isinstance(align, str) and align not in ('center', 'first', 'last'):
        raise ValueError(f'Unknown align method: {align}')
    # Because we change feature locations we have to make a copy of the features.
    # Before, we have to search for the features to be centered.
    pos_first_ft = {}
    stop_pos_last_ft = {}
    pos_align = {}
    for gkey in grouped:
        pos_first_ft[gkey] = min(grouped[gkey], key=lambda ft: ft.loc.start).loc.start
        stop_pos_last_ft[gkey] = max(grouped[gkey], key=lambda ft: ft.loc.stop).loc.stop
        if align == 'first':
            pos_align[gkey] = pos_first_ft[gkey]
        elif align == 'last':
            pos_align[gkey] = stop_pos_last_ft[gkey]
        elif align == 'center':
            pos_align[gkey] = (stop_pos_last_ft[gkey]- pos_first_ft[gkey]) // 2 + pos_first_ft[gkey]
        else:
            for ft in grouped[gkey]:
                if ft in align:
                    pos_align[gkey] = ft.loc.start
                    break
            else:
                warn(f'No feature found for centering group {gkey}')
                pos_align[gkey] = pos_first_ft[gkey]
    if not inplace:
        from copy import deepcopy
        grouped = deepcopy(grouped)
    seq_start = {}
    for gkey in grouped:
        shift = - pos_align[gkey] + max(pos_align.values())
        for ft in grouped[gkey]:
            ft.locs.shift(shift)
        pos_first_ft[gkey] += shift
        stop_pos_last_ft[gkey] += shift
        seq_start[gkey] = shift
    if crop:
        shift = -min(pos_first_ft.values()) + (0 if crop is True else crop)
        for gkey in grouped:
            for ft in grouped[gkey]:
                ft.locs.shift(shift)
            pos_first_ft[gkey] += shift
            stop_pos_last_ft[gkey] += shift
            seq_start[gkey] += shift
    if return_details:
        return grouped, seq_start, pos_first_ft, stop_pos_last_ft
    return grouped


def plot_ftsviewer(
        fts: sugar.FeatureList, fname=None, *,
        seqs=None,
        groupby='seqid',
        align=None, crop=False,
        figsize=None, ncols=1, sharex=False, sharey=False, fig_kw={}, ax=None,
        xticks=False,
        axlabel=True, axlabel_kw={},
        same_colors=True, colorby='type', color=None,
        dpi=None, transparent=None, bbox_inches=None, show=False,
        seqlen=None, first_index=None,
        **kw):
    r"""
    Plot features using DNAFeaturesViewer_ and `~.FeatureList.toftsviewer()`

    :param fts: features
    :param fname: The filename if saving the plot, the figure will be closed afterwards, default: do not save the plot and return figure
    :param seqs: The corresponding sequences, instance of `.BioBasket`,
        will be used to extract the length of each sequence,
        they are also needed when the ``plot_sequence`` option is used.
    :param groupby: Features are grouped by this key, see `~.FeatureList.groupby()`,
        and each group is plotted in a separate ax, defaults to ``'seqid'``.
    :param align: Align features by this method or list of features,
        can be one of the strings ``'center'``, ``'first'``, ``'last'``,
        or a list or BioBasket of features to align each group to the the first found feature in the list.
        By default no alignment is done. If align is used, you should set ``sharex=True``.
    :param crop: If align is used, crop the sequence lines to the first and last feature in each group,
        can be ``True`` to crop to the minimal region including all features,
        or an integer to add some extra space before the first features.
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
    :param seqlen: Can be used to set the sequence length if not using ``seqs``.
    :param first_index: Can be used to set the sequence start index (where the line representing the sequence starts),
        for aligned features this is determined automatically.
        If the parameter is set, it overrides the automatic detection.
        Can be either an int or a dict mapping group keys to first indices.
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
        raise ImportError('Please install dna_features_viewer_lite to use ftsviewer functionality') from ex
    import matplotlib.pyplot as plt
    grouped = fts.groupby(groupby, flatten=True)
    N = len(grouped)
    if N > 1 and same_colors:
        from sugar.imaging.alignment import _get_fts_colordict
        color, colorby = _get_fts_colordict(fts, color, colorby)
    if align:
        if not sharex:
            warn('It is recommended to use sharex=True when align is used')
        grouped, seq_start, pos_first_ft, stop_pos_last_ft = _align_fts(grouped, align, crop=crop, inplace=False, return_details=True)
    else:
        first_index2 = None
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
    xlims = []
    for i, gkey in enumerate(grouped):
        # fi is the position where the sequence aka line begins
        if first_index:
            fi = first_index if isinstance(first_index, int) else first_index[gkey]
        elif align and crop is not False:
            fi = pos_first_ft[gkey] + (0 if crop is True else -crop)
            seqlen = (stop_pos_last_ft[gkey] - pos_first_ft[gkey]) + (0 if crop is True else 2 * crop)
        elif align:
            fi = seq_start[gkey]
        else:
            fi = 0
        if seqs is not None:
            seq = seqs.d[grouped[gkey][0].seqid]
        if axes is not None:
            ax = axes[i % nrows, i // nrows]
        record = grouped[gkey].toftsviewer(seq=seq, colorby=colorby, color=color, seqlen=seqlen,
                                           first_index=fi,
                                           **kw)
        record.plot(ax=ax, **kw2)
        xlims.append(ax.get_xlim())
        if axlabel is True and gkey is not None or axlabel not in (True, False, None):
            axlabel_kw = axlabel_kw.copy()
            axlabel_kw.setdefault('xy', (0, 0.8))
            axlabel_kw.setdefault('xycoords', 'axes fraction')
            axlabel_kw.setdefault('va', 'top')
            l_ = ' '.join(map(str, gkey)) if axlabel is True else grouped[gkey][0].meta[axlabel]
            ax.annotate(l_, **axlabel_kw)
    if sharex:
        x0lim = min(x[0] for x in xlims)
        x1lim = max(x[1] for x in xlims)
        ax.set_xlim(x0lim, x1lim)
    if xticks is True and sharex:
        import matplotlib.ticker as ticker
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    elif xticks and xticks is not True:
        ax.set_xticks(xticks)
    elif not xticks:
        if axes is not None:
            for i in range(ncols*nrows):
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
