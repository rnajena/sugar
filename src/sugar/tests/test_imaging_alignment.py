# (C) 2024, Tom Eulenfeld, MIT license
import pytest
from sugar import read


def test_plot_alignment(outdir):
    """
    Just test that no error occurs
    """
    pytest.importorskip('matplotlib', reason='require matplotlib module')
    import matplotlib.pyplot as plt
    seqs = read().sl(update_fts=True)[:, :100]
    seqs[1][:10] = '-' * 10
    seqs[1][:5] = ' ' * 5
    seqs.fts = seqs.fts.slice(10, 15)[:1]
    fig, axes = plt.subplots(5, figsize=(20, 10))
    seqs.plot_alignment(ax=axes[0], xticks=False)
    seqs.plot_alignment(ax=axes[1], aspect=2, color='flower', symbols=True)
    seqs.plot_alignment(ax=axes[2], fts=True, aspect=2, rasterized=True, color='0.8', symbols=True, symbol_color='flower')
    seqs.plot_alignment(ax=axes[3], fts=True, aspect=2, color='0.8', symbols=True, fts_color='white', fts_alpha=0.5)
    seqs.plot_alignment(ax=axes[4], fts=True, fts_display='box', aspect=2, color='0.8', symbols=True)
    plt.savefig(outdir / 'test_plot_alignment.pdf')
    plt.close()

    fig, ax = plt.subplots(1, figsize=(20, 5))
    seqs.plot_alignment(ax=ax, extent=[0, 100, 0, 10])
    seqs.plot_alignment(ax=ax, extent=[0, 100, -20, -2], rasterized=True, color='flower', symbols=True)
    seqs.plot_alignment(ax=ax, extent=[0, 100, -40, -22], color='0.8', symbols=True, fts=True, fts_alpha=0.5)
    seqs.plot_alignment(ax=ax, extent=[0, 50, -60, -42], fts=True, fts_display='box', aspect=2, color='0.8', symbols=True, fts_color='blue', fts_alpha=0.5)
    plt.savefig(outdir / 'test_plot_alignment2.pdf')
    plt.close()


@pytest.mark.webtest
def test_plot_alignment_examples(outdir):
    pytest.importorskip('matplotlib', reason='require matplotlib module')
    seqs = read('https://osf.io/download/j2wyv')
    seqs.plot_alignment(outdir / 'ali1.png', color='gray', figsize=(10, 2))
    seqs[10:20, 70:120].plot_alignment('ali2.png', color=None, figsize=(10,4),
                                symbols=True, aspect=2, alpha=0.5, xticks=False, bbox_inches='tight')
    seqs2 = seqs[:5, :150].copy()
    seqs2.translate(complete=True).plot_alignment(
        outdir / 'ali3.png', color='flower', figsize=(10,4), symbols=True, aspect=2, alpha=0.5, xticks=False, edgecolors='w')


def _debug_wrong_symbol_size(ax, x, n):
    # error debug message, see issue #6
    fig = ax.get_figure()
    dpiinv = fig.dpi_scale_trans.inverted()
    bbox = ax.get_window_extent().transformed(dpiinv)
    symbol_size = bbox.width * abs((x[-1] - x[0]) / (ax.get_xlim()[1] - ax.get_xlim()[0])) / n * 100
    msg = ('Wrong symbol size\n'
           f'{symbol_size=:.2f} = {bbox.width=:} * abs({(x[-1] - x[0])} / ({ax.get_xlim()[1]=:.2f} - {ax.get_xlim()[0]=:.2f})) / {n} * 100')
    return msg


def test_aspect_symbol_size(outdir):
    mpl = pytest.importorskip('matplotlib', reason='require matplotlib module')
    import matplotlib.pyplot as plt
    seqs = read()[:, :100]

    with mpl.rc_context({'figure.autolayout': False,'figure.constrained_layout.use': False,
                         'figure.subplot.left': 0.1, 'figure.subplot.right':  0.9}):
        for adjustable in ('box', 'datalim'):
            fig = plt.figure(figsize=(12.5, 8), dpi=400)
            ax = fig.add_subplot()
            ax.set_adjustable(adjustable)
            ax.set_xlim(-5, 15)
            ax.set_ylim(-2, 12)
            seqs.plot_alignment(outdir / f'test_plot_ali_aspect1_{adjustable}.png', ax=ax, symbols=True, color=None, aspect=2, extent=[0, 10, 0, 10], show_spines=True, dpi=300)
            assert ax.get_aspect() == 0.04
            text = [child for child in ax.get_children() if hasattr(child, 'get_text') and child.get_text() == 'A'][0]
            assert round(text.get_fontsize(), 1) == 5.0, _debug_wrong_symbol_size(ax, x=(0, 10), n=100)
            plt.close(fig)

        seqs = seqs[:, :2]
        for adjustable in ('box', 'datalim'):
            fig = plt.figure(figsize=(2, 1), dpi=100)
            ax = fig.add_subplot()
            ax.set_adjustable(adjustable)
            ax.set_xlim(-2, 2)
            seqs.plot_alignment(outdir / f'test_plot_ali_aspect2_{adjustable}.png', ax=ax, symbols=True, color=None, aspect=2, show_spines=True, scale_symbol_size=1/1.925, dpi=100)
            assert ax.get_aspect() == 2
            text = [child for child in ax.get_children() if hasattr(child, 'get_text') and child.get_text() == 'A'][0]
            assert round(text.get_fontsize(), 1) == 10.0, _debug_wrong_symbol_size(ax, x=(-0.5, 1.5), n=2*1.925)
            plt.close(fig)
