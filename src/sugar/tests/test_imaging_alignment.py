# (C) 2024, Tom Eulenfeld, MIT license
import pytest
from sugar import read
from sugar.tests.util import _changetmpdir


def test_plot_alignment():
    """
    Just test that no error occurs, you can check the image by uncommenting the last line
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
    # plt.savefig('test_plot_alignment.pdf')
    plt.close()

    fig, ax = plt.subplots(1, figsize=(20, 5))
    seqs.plot_alignment(ax=ax, extent=[0, 100, 0, 10])
    seqs.plot_alignment(ax=ax, extent=[0, 100, -20, -2], rasterized=True, color='flower', symbols=True)
    seqs.plot_alignment(ax=ax, extent=[0, 100, -40, -22], color='0.8', symbols=True, fts=True, fts_alpha=0.5)
    seqs.plot_alignment(ax=ax, extent=[0, 50, -60, -42], fts=True, fts_display='box', aspect=2, color='0.8', symbols=True, fts_color='blue', fts_alpha=0.5)
    # plt.savefig('test_plot_alignment2.pdf')
    plt.close()


@pytest.mark.webtest
def test_plot_alignment_examples():
    pytest.importorskip('matplotlib', reason='require matplotlib module')
    seqs = read('https://osf.io/download/j2wyv')
    with _changetmpdir() as tmpdir:
        tmpdir = './'
        seqs.plot_alignment(tmpdir + 'ali1.png', figsize=(10, 2))
        seqs[10:20, 70:120].plot_alignment(tmpdir + 'ali2.png', color=None, figsize=(10,4),
                                    symbols=True, aspect=2, alpha=0.5, xticks=False, bbox_inches='tight')
        seqs2 = seqs[:5, :150].copy()
        seqs2.translate(complete=True).plot_alignment(
            tmpdir + 'ali3.png', color='flower', figsize=(10,4), symbols=True, aspect=2, alpha=0.5, xticks=False, edgecolors='w')




