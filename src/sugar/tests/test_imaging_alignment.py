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


def test_aspect_symbol_size(outdir):
    pytest.importorskip('matplotlib', reason='require matplotlib module')
    import matplotlib.pyplot as plt
    seqs = read()[:, :100]

    for adjustable in ('box', 'datalim'):
        fig = plt.figure(figsize=(10, 8), dpi=100)
        ax = fig.add_subplot()
        ax.set_adjustable(adjustable)
        ax.set_xlim(-10, 150)
        ax.set_xlim(-2, 15)
        seqs.plot_alignment(outdir / f'test_plot_ali_aspect1_{adjustable}.png', ax=ax, symbols=True, color=None, aspect=2, extent=[0, 10, 0, 10], show_spines=True)
        assert ax.get_aspect() == 0.04
        text = [child for child in ax.get_children() if hasattr(child, 'get_text') and child.get_text() == 'A'][0]
        assert round(text.get_fontsize(), 1) == 4.6
        plt.close(fig)

    seqs = seqs[:, :2]
    for adjustable in ('box', 'datalim'):
        fig = plt.figure(figsize=(2, 1), dpi=100)
        ax = fig.add_subplot()
        ax.set_adjustable(adjustable)
        ax.set_xlim(-1, 5)
        seqs.plot_alignment(outdir / f'test_plot_ali_aspect2_{adjustable}.png', ax=ax, symbols=True, color=None, aspect=2, show_spines=True)
        assert ax.get_aspect() == 2
        text = [child for child in ax.get_children() if hasattr(child, 'get_text') and child.get_text() == 'A'][0]
        assert round(text.get_fontsize(), 2) == 19.25
        plt.close(fig)
