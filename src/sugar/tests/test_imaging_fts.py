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
    seq = client.get_seq('AF086833')
    fts = seq.fts.select('cds')
    for ft in fts:
        ft.meta.name = ft.meta._genbank.gene
    fts.plot_ftsviewer(outdir / 'fts3.png', colorby='name', seqlen=len(seq), figsize=(6, 2.5))
    orfs = seq.find_orfs(len_ge=500)
    orfs.plot_ftsviewer(outdir / 'fts4.png', label=None, colorby='rf', seqlen=len(seq), figsize=(6, 2.5))

    seq2 = seq.copy()
    seq2.data = seq2.data[:15_000]
    seq2.id = 'fake'
    fts2 = fts.copy()[:-1]
    for ft in fts2:
        ft.seqid = 'fake'
        ft.locs.shift(1_000)
    fts3 = fts + fts2
    align = [fts[2], fts2[2]]
    kw = dict(colorby='name', figsize=(6, 2.5), sharex=True, xticks=True)
    fts3.plot_ftsviewer(outdir / 'fts_align.png', align=align, **kw)
    from sugar import BioBasket
    seqs3 = BioBasket([seq, seq2])
    seqs3.fts = fts3
    seqs3.plot_ftsviewer(outdir / 'fts_align_seqs_center.png', align='center', **kw)
    seqs3.plot_ftsviewer(outdir / 'fts_align_seqs_nothing_done.png',  **kw)
    seqs3.plot_ftsviewer(outdir / 'fts_align_seqs.png', align=align, **kw)
    seqs3.plot_ftsviewer(outdir / 'fts_align_seqs_first.png', align='first', **kw)
    seqs3.plot_ftsviewer(outdir / 'fts_align_seqs_center_crop100.png', align='center', crop=100, **kw)
    seqs3.plot_ftsviewer(outdir / 'fts_align_seqs_last_crop.png', align='last', crop=True, **kw)
    # test hatch and alternative way to set color, hatch, label
    seqs3.fts[-1].meta._ftsviewer_hatch='.'
    seqs3.fts[-1].meta._ftsviewer_color='C0'
    seqs3.fts[-1].meta._ftsviewer_label='test'
    try:
        seqs3.plot_ftsviewer(outdir / 'fts_align_seqs_crop.png', align=align, crop=True, **kw)
    except AttributeError:
        pytest.xfail('old dna_features_viewer_lite version or original dna_features_viewer package without hatch support')
