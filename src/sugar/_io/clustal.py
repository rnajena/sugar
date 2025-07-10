# (C) 2024, Tom Eulenfeld, MIT license
"""
`Clustal`_ IO
"""
from warnings import warn
from sugar import BioBasket, BioSeq, __version__
from sugar._io.util import _add_fmt_doc


def is_clustal(f, **kw):
    content = f.read(7)
    return content == 'CLUSTAL'


@_add_fmt_doc('read')
def read_clustal(f):
    """
    Read Clustal alignment from file into `.BioBasket`
    """
    seqs = {}
    for line in f:
        line = line.strip()
        if line == '' or line.startswith('CLUSTAL'):
            continue
        parts = line.split()
        if len(parts) >= 2:
            id_, seqdata = parts[:2]
            seqs.setdefault(id_, []).append(seqdata)
        else:
            # do not read conservation line
            # if we do this in the future we cannot rely on split(),
            # because characters at the beginning or end of the conservation string might be whitespaces
            pass
    return BioBasket([BioSeq(''.join(data), id=id_) for id_, data in seqs.items()])


_SEP = 5
_CHARS = 60


@_add_fmt_doc('write')
def write_clustal(seqs, f, header_sugar=True, header=None):
    """
    Write `.BioBasket` to Clustal format

    :param bool header_sugar: Append header with sugar version to the first line, default is True
    :param str header: More information appended to the first line
    """
    content = []
    line = 'CLUSTAL'
    if header_sugar:
        line = line + f' format written by sugar v{__version__}'
    if header:
        line = line + ' ' + header
    content = [line + '\n\n']
    idlen = max(len(id_) for id_ in seqs.ids)
    lens = {len(seq) for seq in seqs}
    if len(lens) > 1:
        warn('Writing Clustal file with sequences of different lengths')
    if len(seqs) > 0:
        pos = 0
        while pos < max(lens):
            for seq in seqs:
                content.append(f'{seq.id:<{idlen+_SEP}}{seq.data[pos:pos+_CHARS]}\n')
            content.append('\n')
            pos += _CHARS
    f.write(''.join(content))
