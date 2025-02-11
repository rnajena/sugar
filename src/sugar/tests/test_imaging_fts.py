# (C) 2024, Tom Eulenfeld, MIT license
import pytest
from sugar import read
from sugar.tests.util import _changetmpdir

# set path to save images
_PATH = None


def test_plot_fts_via_ftsviewer():
    pytest.importorskip('matplotlib', reason='require matplotlib module')
    pytest.importorskip('dna_features_viewer', reason='require dna_features_viewer module')
    import matplotlib.pyplot as plt
    seqs = read()
    record = seqs[1].toftsviewer()
    record.plot()
    with _changetmpdir(_PATH) as tmpdir:
        plt.savefig('fts1.png')
    plt.close()


@pytest.mark.webtest
def test_plot_fts_via_ftsviewer_genes():
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
    record = fts.toftsviewer(colorby='name', seqlen=len(seq))
    record.plot()
    with _changetmpdir(_PATH) as tmpdir:
        plt.savefig('fts2.png')
    plt.close()
