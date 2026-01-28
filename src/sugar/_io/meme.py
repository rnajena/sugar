# (C) 2024, Tom Eulenfeld, MIT license
"""
Read support for file formats associated with the `MEME`_ Suite
"""


from sugar._io.util import _add_fmt_doc
from sugar.core.fts import Feature, FeatureList


def is_fts_meme_txt(f, **kw):
    content = f.read(150)
    lines = content.splitlines()
    return '*' in lines[0] and lines[1].startswith('MEME')


def _read_motifs_meme_txt(f):
    """
    Parse MEME txt file and return dictionary of motifs
    """
    content = f.read()
    lines = content.splitlines()
    motif = None
    motifs = []
    important_part = False
    for line in lines:
        if line.startswith('MOTIF'):
            if motif:
                motifs.append(motif)
            _, motif_seq, motif_id, *_, evalue = line.split()
            motif = {'motif': motif_seq, 'motif_id': motif_id, 'evalue': float(evalue), 'hits': []}
        elif motif:
            if important_part:
                if line.strip() == '':
                    important_part = False
                elif not line.startswith('---') and not line.startswith('Sequence name'):
                    seqid, strand, start, pvalue, *_ = line.split()
                    motif['hits'].append({'seqid': seqid, 'strand': strand, 'start': int(start)-1, 'pvalue': float(pvalue)})
            elif 'sites sorted by position p-value' in line:
                important_part = True
                assert len(motif['hits']) == 0
    if motif:
        motifs.append(motif)
    return motifs


#https://regex101.com/r/82v1Zy/2
_MEME_TXT_REGEX = r"\*+\nMOTIF\s+(?P<motif>\w+)\s+(?P<motif_no>[\w-]+).+E-value = (?P<evalue>[\d\.e+-]+)\s*\n\*+\n.*Motif\s+(?P=motif)\s+(?P=motif_no) sites sorted by position p-value\n-+\nSeq[\s\w-]+\n[-\s]+\n(?P<hitdata>.*?)\n-+\n"

def _read_motifs_meme_txt_regex(f):
    """
    Parse MEME txt file with regex and return dictionary of motifs
    """
    import re
    content = f.read()
    motifs = []
    for match in re.finditer(_MEME_TXT_REGEX, content, re.DOTALL):
        motif = {
            'motif': match.group('motif'),
            'motif_id': match.group('motif_no'),
            'evalue': float(match.group('evalue')),
            'hits': []
        }
        hitdata = match.group('hitdata')
        for line in hitdata.splitlines():
            seqid, strand, start, pvalue, *_ = line.split()
            motif['hits'].append({'seqid': seqid, 'strand': strand, 'start': int(start)-1, 'pvalue': float(pvalue)})
        motifs.append(motif)
    return motifs


@_add_fmt_doc('read_fts')
def read_fts_meme_txt(f, *, ftype=None, engine='regex'):
    """
    MEME txt reader

    :param str ftype: Feature type of returned features
    :param str engine: Select ``'txt'`` or ``'regex'`` based parsing method, both should give the same result,
        default is ``'regex'``

    Returns a flat `.FeatureList` of all motif hits.
    To get a dictionary of motifs call `FeatureList.groupby('motif')<.FeatureList.groupby>` or
    `FeatureList.groupby('motif_id')<.FeatureList.groupby>` on the returned object.
    Alternatively, use the `_read_motifs_meme_txt_regex` function directly.
    """
    if engine == 'regex':
        motifs = _read_motifs_meme_txt_regex(f)
    elif engine == 'txt':
        motifs = _read_motifs_meme_txt(f)
    else:
        raise ValueError("Unknown read engine, must be one of 'regex', 'txt'")
    fts = []
    for motif in motifs:
        for hit in motif['hits']:
            ft = Feature(ftype,
                         start=hit['start'], stop=hit['start'] + len(motif['motif']),
                         strand=hit['strand'],
                         meta={
                             'seqid': hit['seqid'],
                             'motif': motif['motif'],
                             'motif_id': motif['motif_id'],
                             'evalue': motif['evalue'],
                             '_meme': {'pvalue': hit['pvalue']}
                         })
            fts.append(ft)
    return FeatureList(fts)
