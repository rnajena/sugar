# (C) 2024, Tom Eulenfeld, MIT license
"""
FASTA IO
"""
import re

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


CHS = '[^,|;\s]'

# capture some cases from https://de.wikipedia.org/wiki/FASTA-Format
IDPATTERN= (
    f'[^\s]*gb[:|]({CHS}+)'  # gb:id, gb|id
    f'|[^\s]*(?:emb|dbj|sp|tr|ref|lcl)[|]({CHS}+)' # xxx|id
    f'|({CHS}+)'  # "id ", "id;", "id|", "id,"
    )


def _id_from_header(header):
    id_ = None
    if header != '':
        match = re.match(IDPATTERN, header)
        if match is not None:
            for id_ in match.groups():
                if id_ is not None:
                    break
    return id_


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
            id_ = _id_from_header(header)
            # if line == '':
            #     id_ = None
            # elif ' ' not in line and '|' not in line:
            #     id_ = line
            # elif '|' in line:
            #     id_ = line.split('|')[-2]
            # else:
            #     id_ = line.split(maxsplit=1)[0]
            data = []
        elif line.startswith(';'):  # line is a comment
            pass
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
        header = seq.meta._fasta.header.removeprefix(id_)
        if header.strip() == '':
            header = ''
        else:
            header = ' ' + header
    else:
        header = ''
    content = f'>{id_}{header}\n{seq.data}\n'
    f.write(content)
