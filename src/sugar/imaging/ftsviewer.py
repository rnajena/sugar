import sugar
from sugar.core.util import _pop_kws_for_func


def plot_ftsviewer(
        fts: sugar.FeatureList, fname=None, *,
        seqs=None,
        groupby='seqid', axlabel=True,
        figsize=None, fig_kw={}, axlabel_kw={},
        ncols=1,
        dpi=None, transparent=None, bbox_inches=None, show=False,
        **kw):
    import matplotlib.pyplot as plt
    from dna_features_viewer import GraphicRecord
    grouped = fts.groupby(groupby, flatten=True)
    N = len(grouped)
    nrows = (N + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, squeeze=False, figsize=figsize, **fig_kw)
    seq = None
    kw2 = _pop_kws_for_func(kw, GraphicRecord.plot)
    for i, gkey in enumerate(grouped):
        ax = axes[i % nrows, i // nrows]
        if seqs is not None:
            seq = seqs.d[grouped[gkey][0].seqid]
        record = grouped[gkey].toftsviewer(seq=seq, **kw)
        record.plot(ax=ax, **kw2)
        if gkey is not None:
            axlabel_kw = axlabel_kw.copy()
            axlabel_kw.setdefault('xy', (0, 0.5))
            axlabel_kw.setdefault('xycoords', 'axes fraction')
            ax.annotate(' '.join(map(str, gkey)), **axlabel_kw)
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
