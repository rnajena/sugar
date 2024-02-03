# (C) 2024, Tom Eulenfeld, MIT license

from sugar import BioSeq


EXT = ['fasta']

def is_format(f):
    content = f.read(50)
    return content.strip().startswith('>')


def create_bioseq(datalines, id_, header):
    data = ''.join(datalines)
    seq = BioSeq(data, id=id_)
    seq.meta._fasta = {}
    seq.meta._fasta.header = header
    return seq


def iter_(f):
    content = f.read()
    id_ = None
    header = None
    data = None
    for line in content.splitlines():
        if line.startswith('>'):
            if data is not None:
                yield create_bioseq(data, id_, header)
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
        yield create_bioseq(data, id_, header)


def append(seq, f, wrap=None):
    id_ = seq.id or ''
    if '_fasta' in seq.meta and 'comment' in seq.meta._fasta:
        comment = ' ' + seq.meta._fasta.comment
    else:
        comment = ''
    content = f'> {id_}{comment}\n{seq.data}\n'
    if wrap:
        from textwrap import fill
        content = fill(content, width=wrap)
    f.write(content)
