# (C) 2024, Tom Eulenfeld, MIT license
"""
`Stockholm`_ IO
"""

from warnings import warn

from sugar.core.seq import Attr, BioSeq, BioBasket
from sugar._io.util import _add_fmt_doc


filename_extensions_stockholm = ['stk', 'sto', 'stockholm']


def is_stockholm(f, **kw):
    content = f.read(11)
    return content == '# STOCKHOLM'


def row2fts(row, type_=None, seqid=None, fillchars='-.', splitchar='|'):
    import re
    from sugar.core.fts import Defect, Feature, FeatureList, Location
    i = 0
    fts = []
    while i < len(row):
        j = row.find(splitchar, i+1)
        if j == -1:
            j = len(row)
        names = set(re.split(f'[{fillchars}]+', row[i+1:j])) - {''}
        if len(names) > 0:
            stop = min(j+1, len(row))
            if len(names) > 1:
                warn('More than one name for feature at location {i}:{stop}, use shortest name')
                names = sorted(names, key=lambda n: len(n))[:1]
            name = names.pop()
            if i == 0 and row[0] != splitchar and j == len(row):
                assert row[-1] != splitchar
                defect = Defect.MISS_LEFT | Defect.MISS_RIGHT
            if i == 0 and row[0] != splitchar:
                defect = Defect.MISS_LEFT
            elif j == len(row):
                assert row[-1] != splitchar
                defect = Defect.MISS_RIGHT
            else:
                defect = Defect.NONE
            loc = Location(i, stop, defect=defect)
            ft = Feature(type_, [loc])
            ft.meta.name = name
            if seqid is not None:
                ft.seqid = seqid
            fts.append(ft)
        i = j
    return FeatureList(fts)


def fts2row(fts, inchar='-', outchar='.', splitchar='|', lensec=None):
    openchar = closechar = bothchar = splitchar
    row = []
    last_stop = 0
    for ft in fts.sort():
        if len(ft.locs) > 1:
            warn('More than one location in feature, use full loc_range')
        start, stop = ft.locs.range
        dif = start - last_stop
        if dif > 0:
            row.append(outchar * dif)
        elif dif == -1:  # use the same | for last and this feature
            assert row[-1][-1] == closechar
            row[-1] = row[-1][:-1]
        elif dif < -1:
            raise ValueError('Features overlap more than one residue')
        l = stop - start
        name = '' if ft.name is None else ft.name
        rowp = name
        while l > len(name) + len(rowp) + 150:
            rowp = name + inchar * 100 + rowp
        rowp = rowp.center(l, inchar)
        if l < len(name)-2:
            warn('Feature name too long, shorten it')
            rowp = inchar + name[:l-2] + inchar
        if (ft.loc.defect.MISS_LEFT not in ft.loc.defect and
                ft.loc.defect.BEYOND_LEFT not in ft.loc.defect):
            rowp = (bothchar if dif == -1 else openchar) + rowp[1:]
        if (ft.loc.defect.MISS_RIGHT not in ft.locs[-1].defect and
                ft.loc.defect.BEYOND_RIGHT not in ft.locs[-1].defect):
            rowp = rowp[:-1] + closechar
        assert len(rowp) == l
        row.append(rowp)
        last_stop = stop
    return ''.join(row).ljust(lensec or len(row), outchar)


@_add_fmt_doc('read')
def read_stockholm(f, comments=None):
    """
    Read Stockholm file

    :param list comments: comment lines inside the file are stored in
        the comments list (optional)

    .. note::
        By default only the first alignment is read and returned.
        If your file contains more than one alignment, you can read
        some or all of them with::

            from sugar import read
            with open('example_stockholm_multi.stk') as f:
                seqs1 = read(f, 'stockholm')  # read 1st alignment
                seqs2 = read(f, 'stockholm')  # read 2nd alignment
    """
    seqs = []
    gf = Attr()
    gc = Attr()
    gs = {}
    gr = {}
    seqs = {}
    for line in f:
        line = line.strip()
        if line == '' or line.startswith('# STOCKHOLM'):
            continue
        elif line.startswith('#=GF CC written by sugar'):
            continue
        elif line.startswith('#=GF'):
            _, key, val = line.split(maxsplit=2)
            gf[key] = gf[key] + ' ' + val if key in gf else val
        elif line.startswith('#=GC'):
            _, key, val = line.split(maxsplit=2)
            gc[key] = gc.get(key, '') + val
        elif line.startswith('#=GS'):
            _, seqid, key, val = line.split(maxsplit=3)
            gs.setdefault(seqid, {})
            gs[seqid][key] = gs[seqid][key] + ' ' + val if key in gs[seqid] else val
        elif line.startswith('#=GR'):
            _, seqid, key, val = line.split(maxsplit=3)
            gr[seqid][key] = gr.setdefault(seqid, {}).get(key, '') + val
        elif line.startswith('#'):  # ignore other comments
            if comments is not None:
                comments.append(line)
        elif line.startswith('//'):  # end of alignment, stop reading
            break
        elif ' ' in line:  # sequence
            key, val = line.split(maxsplit=1)
            seqs[key] = seqs.get(key, '') + val
            # if seq is not None and seq.id == key:
            #     seq += val
            # else:
            #     if seq is not None:
            #         seqs.append(seq)
            #     seq = BioSeq(val, id=key)
        else:  # should not happen
            raise ValueError('Invalid stockholm format, please contact devs')
    seqs = [BioSeq(val, id=key) for key, val in seqs.items()]
    for seq in seqs:
        if seq.id in gs:
            seq.meta.setdefault('_stockholm', Attr()).GS = gs[seq.id]
        if seq.id in gr:
            seq.meta.setdefault('_stockholm', Attr()).GR = gr[seq.id]
    seqs = BioBasket(seqs)
    if len(gf) > 0:
        seqs.meta.setdefault('_stockholm', Attr()).GF = gf
    if len(gc) > 0:
        seqs.meta.setdefault('_stockholm', Attr()).GC = gc
    return seqs


@_add_fmt_doc('write')
def write_stockholm(seqs, f, header_sugar=True):
    """
    Write sequences to Stockholm file

    :param bool header_sugar: Add a comment to the header with sugar version, default is True
    """
    lines = ['# STOCKHOLM 1.0']
    if header_sugar:
        from sugar import __version__
        lines.append(f'#=GF CC written by sugar v{__version__}')
    try:
        gf = seqs.meta._stockholm.GF
    except AttributeError:
        pass
    else:
        for key, value in gf.items():
            lines.append(f'#=GF {key} {value}')
    for seq in seqs:
        try:
            gs = seq.meta._stockholm.GS
        except (AttributeError, KeyError):
            pass
        else:
            for key, value in gs.items():
                lines.append(f'#=GS {seq.id} {key} {value}')
    for seq in seqs:
        lines.append(f'{seq.id} {seq}')
        try:
            gr = seq.meta._stockholm.GR
        except (AttributeError, KeyError):
            pass
        else:
            for key, value in gr.items():
                lines.append(f'#=GR {seq.id} {key} {value}')
    try:
        gc = seqs.meta._stockholm.GC
    except AttributeError:
        pass
    else:
        for key, value in gc.items():
            lines.append(f'#=GC {key} {value}')
    lines.append('//\n')
    f.write('\n'.join(lines))
