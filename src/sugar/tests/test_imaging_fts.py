# (C) 2024, Tom Eulenfeld, MIT license
import pytest
from sugar import read
from contextlib import chdir


def test_plot_ftsviewer(outdir):
    pytest.importorskip('matplotlib', reason='require matplotlib module')
    pytest.importorskip('dna_features_viewer', reason='require dna_features_viewer module')
    import matplotlib.pyplot as plt
    seqs = read()
    seqs.plot_ftsviewer(outdir / 'fts1.png')
    fig = seqs[0].plot_ftsviewer()
    plt.close(fig)
    ax = plt.axes()
    fig = seqs[0].plot_ftsviewer(ax=ax)
    plt.close(fig)


def test_plot_ftsviewer_record(outdir):
    pytest.importorskip('matplotlib', reason='require matplotlib module')
    pytest.importorskip('dna_features_viewer', reason='require dna_features_viewer module')
    import matplotlib.pyplot as plt
    seqs = read()
    record = seqs[1].toftsviewer()
    record.plot()
    plt.savefig(outdir / 'fts2.png')
    plt.close()


@pytest.mark.webtest
def test_plot_ftsviewer_genes(outdir):
    pytest.importorskip('matplotlib', reason='require matplotlib module')
    pytest.importorskip('dna_features_viewer', reason='require dna_features_viewer module')
    import matplotlib.pyplot as plt
    from sugar.web import Entrez
    client = Entrez()
    with pytest.warns():  # join in genbank loc is untested
        seq = client.get_seq('AF086833')
    fts = seq.fts.select('cds')
    for ft in fts:
        ft.meta.name = ft.meta._genbank.gene
    fts.plot_ftsviewer(outdir / 'fts3.png', colorby='name', seqlen=len(seq), figsize=(6, 2.5))
    orfs = seq.find_orfs(len_ge=500)
    orfs.plot_ftsviewer(outdir / 'fts4.png', label=None, colorby='rf', seqlen=len(seq), figsize=(6, 2.5))
    plt.close()
