# (C) 2024, Tom Eulenfeld, MIT license
import pytest
from sugar import read
from sugar.tests.util import _changetmpdir

# set path to save images
_PATH = None


def test_plot_ftsviewer():
    pytest.importorskip('matplotlib', reason='require matplotlib module')
    pytest.importorskip('dna_features_viewer', reason='require dna_features_viewer module')
    import matplotlib.pyplot as plt
    seqs = read()
    with _changetmpdir(_PATH) as tmpdir:
        seqs.plot_ftsviewer('fts1.png')
    fig = seqs[0].plot_ftsviewer()
    plt.close(fig)
    ax = plt.axes()
    fig = seqs[0].plot_ftsviewer(ax=ax)
    plt.close(fig)


def test_plot_ftsviewer_record():
    pytest.importorskip('matplotlib', reason='require matplotlib module')
    pytest.importorskip('dna_features_viewer', reason='require dna_features_viewer module')
    import matplotlib.pyplot as plt
    seqs = read()
    record = seqs[1].toftsviewer()
    record.plot()
    with _changetmpdir(_PATH) as tmpdir:
        plt.savefig('fts2.png')
    plt.close()


@pytest.mark.webtest
def test_plot_ftsviewer_genes():
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
    with _changetmpdir(_PATH) as tmpdir:
        fts.plot_ftsviewer('fts3.png', colorby='name', seqlen=len(seq), figsize=(6, 2.5))
        orfs = seq.find_orfs(rf='all').select(len_gt=500)
        orfs.plot_ftsviewer('fts4.png', label=None, colorby='rf', seqlen=len(seq), figsize=(6, 2.5))
    plt.close()
