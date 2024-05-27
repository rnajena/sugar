# (C) 2024, Tom Eulenfeld, MIT license
"""
FASTA IO
"""

from sugar import BioSeq
from sugar._io.util import _add_fmt_doc


filename_extensions = ['fasta', 'fa']

def is_format(f, **kw):
    content = f.read(50)
    return content.strip().startswith('>')


def _create_bioseq(datalines, id_, header):
    data = ''.join(datalines)
    seq = BioSeq(data, id=id_)
    seq.meta._fasta = {}
    seq.meta._fasta.header = header
    return seq


@_add_fmt_doc('read')
def iter_(f):
    """
    Iterate through a FASTA file and yield `.BioSeq` sequences
    """
    id_ = None
    header = None
    data = None
    for line in f:
        if line.startswith('>'):
            if data is not None:
                yield _create_bioseq(data, id_, header)
            line = line.lstrip('>').strip()
            header = line
            if line == '':
                id_ = None
            elif ' ' not in line and '|' not in line:
                id_ = line
            elif '|' in line:
                id_ = line.split('|')[-2]
            else:
                id_ = line.split(maxsplit=1)[0]
            data = []
        else:
            data.append(line.strip())
    if data is not None:
        yield _create_bioseq(data, id_, header)


@_add_fmt_doc('write')
def append(seq, f):
    """
    Append a `.BioSeq` sequence to a FASTA file
    """
    id_ = seq.id or ''
    if '_fasta' in seq.meta and 'header' in seq.meta._fasta:
        header = ' ' + seq.meta._fasta.header
    else:
        header = ''
    content = f'>{id_}{header}\n{seq.data}\n'
    f.write(content)
