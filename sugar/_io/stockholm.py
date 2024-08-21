# (C) 2024, Tom Eulenfeld, MIT license
# https://en.wikipedia.org/wiki/Stockholm_format
# TODO: convert feature GC line to FeatureList and vice versa
"""
Stockholm IO
"""

from warnings import warn

from sugar.core.seq import Attr, BioSeq, BioBasket
from sugar._io.util import _add_fmt_doc


filename_extensions = ['stk', 'stockholm']


def is_format(f, **kw):
    content = f.read(11)
    return content == '# STOCKHOLM'


def row2fts(row, type_=None, seqid=None):
    import re
    from sugar.core.fts import Feature, FeatureList, Location
    i = 0
    fts = []
    while i < len(row):
        j = row.find('|', i+1)
        if j == -1:
            j = len(row)
        names = set(re.split('[.]+', row[i+1:j])) - {''}
        if len(names) > 0:
            stop = min(j+1, len(row))
            if len(names) > 1:
                warn('More than one name for feature at location {i}:{stop}, use shortest name')
                names = sorted(names, key=lambda n: len(n))[:1]
            name = names.pop()
            if i == 0 and row[0] != '|' and j == len(row):
                assert row[-1] != '|'
                defect = Location.Defect.MISS_LEFT | Location.Defect.MISS_RIGHT
            if i == 0 and row[0] != '|':
                defect = Location.Defect.MISS_LEFT
            elif j == len(row):
                assert row[-1] != '|'
                defect = Location.Defect.MISS_RIGHT
            else:
                defect = Location.Defect.NONE
            loc = Location(i, stop, defect=defect)
            ft = Feature(type_, [loc])
            ft.meta.name = name
            if seqid is not None:
                ft.seqid = seqid
            fts.append(ft)
        i = j
    return FeatureList(fts)


def fts2row(fts):
    row = []
    last_stop = None
    for ft in fts.sort():
        if len(ft.locs) > 1:
            warn('More than one location in feature, use full loc_range')
        start, stop = ft.loc_range
        if last_stop is not None:
            dif = start - last_stop
            if dif > 0:
                row.append('.' * dif)
            elif dif == -1:  # use the same | for last and this feature
                assert row[-1][-1] == '|'
                row[-1] = row[-1][:-1]
            elif dif < -1:
                raise ValueError('Features overlap more than one residue')
        l = stop - start
        name = '' if ft.name is None else ft.name
        rowp = name
        while l > len(name) + len(rowp) + 150:
            rowp = name + '.' * 100 + rowp
        rowp = rowp.center(l, '.')
        if l < len(name)-2:
            warn('Feature name too long, shorten it')
            rowp = '.' + name[:l-2] + '.'
        if ft.loc.Defect.MISS_LEFT not in ft.loc.defect:
            rowp = '|' + rowp[1:]
        if ft.loc.Defect.MISS_RIGHT not in ft.locs[-1].defect:
            rowp = rowp[:-1] + '|'
        assert len(rowp) == l
        row.append(rowp)
        last_stop = stop
    return ''.join(row)


@_add_fmt_doc('read')
def read(f):
    """
    Read Stockholm file
    """
    seqs = []
    gf = Attr()
    gc = Attr()
    gs = {}
    gr = {}
    seq = None
    for line in f:
        line = line.strip()
        if line == '' or line.startswith('# STOCKHOLM'):
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
            continue
        elif line.startswith('//'):  # end of alignment, stop reading
            break
        elif ' ' in line:  # sequence
            key, val = line.split(maxsplit=1)
            if seq is not None and seq.id == key:
                seq += val
            else:
                if seq is not None:
                    seqs.append(seq)
                seq = BioSeq(val, id=key)
        else:  # should not happen
            raise ValueError('Invalid stockholm format, please contact devs')
    if seq is not None:
        seqs.append(seq)
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
def write(seqs, f):
    """
    Write sequences to Stockholm file
    """
    lines = ['# STOCKHOLM 1.0']
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
