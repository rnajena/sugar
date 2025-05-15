# (C) 2024, Tom Eulenfeld, MIT license
"""
`FASTA`_ IO
"""
import re

from sugar import BioSeq
from sugar._io.util import _add_fmt_doc


filename_extensions_fasta = ['fasta', 'fa', 'fna', 'faa', 'fas']

def is_fasta(f, **kw):
    content = f.read(50)
    return content.strip().startswith('>')


def _create_bioseq(datalines, id_, header):
    data = ''.join(datalines)
    seq = BioSeq(data, id=id_)
    seq.meta._fasta = {}
    seq.meta._fasta.header = header
    return seq


CHS = r'[^,|;\s]'  # allowed characters in the seq ID

# capture some cases from https://de.wikipedia.org/wiki/FASTA-Format
IDPATTERN= (
    rf'[^\s]*gb[:|]({CHS}+)'  # gb:id, gb|id
    rf'|[^\s]*(?:emb|dbj|sp|tr|ref|lcl)[|]({CHS}+)' # xxx|id
    rf'|({CHS}+)'  # "id ", "id;", "id|", "id,"
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
def iter_fasta(f, comments=None):
    """
    Iterate through a FASTA file and yield `.BioSeq` sequences

    :param list comments: comment lines inside the file are stored in
        the comments list (optional)
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
            if comments is not None:
                comments.append(line)
        else:
            data.append(line.strip())
    if data is not None:
        yield _create_bioseq(data, id_, header)


@_add_fmt_doc('write')
def append_fasta(seq, f):
    """
    Append a `.BioSeq` sequence to a FASTA file
    """
    id_ = seq.id or ''
    if '_fasta' in seq.meta and 'header' in seq.meta._fasta:
        header = (' ' + seq.meta._fasta.header.removeprefix(id_).lstrip()).rstrip()
    else:
        header = ''
    content = f'>{id_}{header}\n{seq.data}\n'
    f.write(content)
