# (C) 2024, Tom Eulenfeld, MIT license
# https://en.wikipedia.org/wiki/Stockholm_format
# TODO: convert feature GC line to FeatureList and vice versa

from sugar.core.seq import Attr, BioSeq, BioBasket


EXT = ['stk', 'stockholm']


def is_format(f):
    content = f.read(11)
    return content == '# STOCKHOLM'


def read(f):
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
            seq.meta.setdefault('stk', Attr()).GS = gs[seq.id]
        if seq.id in gr:
            seq.meta.setdefault('stk', Attr()).GR = gr[seq.id]
    seqs = BioBasket(seqs)
    if len(gf) > 0:
        seqs.meta.setdefault('stk', Attr()).GF = gf
    if len(gc) > 0:
        seqs.meta.setdefault('stk', Attr()).GC = gc
    return seqs


def write(seqs, f):
    lines = ['# STOCKHOLM 1.0']
    try:
        gf = seqs.meta.stk.GF
    except AttributeError:
        pass
    else:
        for key, value in gf.items():
            lines.append(f'#=GF {key} {value}')
    for seq in seqs:
        try:
            gs = seq.meta.stk.GS
        except (AttributeError, KeyError):
            pass
        else:
            for key, value in gs.items():
                lines.append(f'#=GS {seq.id} {key} {value}')
    for seq in seqs:
        lines.append(f'{seq.id} {seq}')
        try:
            gr = seq.meta.stk.GR
        except (AttributeError, KeyError):
            pass
        else:
            for key, value in gr.items():
                lines.append(f'#=GR {seq.id} {key} {value}')
    try:
        gc = seqs.meta.stk.GC
    except AttributeError:
        pass
    else:
        for key, value in gc.items():
            lines.append(f'#=GC {key} {value}')
    lines.append('//\n')
    f.write('\n'.join(lines))
