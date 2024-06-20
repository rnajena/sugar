# (C) 2024, Tom Eulenfeld, MIT license

from sugar.data import gcode
import warnings


def translate(seq, *, complete=False, check_start=None, check_stop=False,
              warn=False, astop='X', gap='-', gap_after=2, tt=1):
    """
    Translate a string or `.BioSeq` object into an amino acid string

    :param bool complete: If set to ``True`` ignores stop codons,
        otherwise the translation is stopped before the first stop codon
    :param bool check_start: Check that the first codon is a start codon,
        default is False for ``complete=False`` otherwise True
    :param bool check_stop: Check that the sequence ends with the first stop
        codon, default is False
    :param bool warn: Warn if the first codon might not be a start codon
        (if ``check_start=True``) and warn for amigious stop codons,
        default is False
    :param str astop: Symbol for ambigious stop codons
    :param str gap: gap character, default ``'-'``, set to ``None``
       to raise an error for non nucleotide characters
    :param int gap_after: A single gap in the amino acis string is
        written after the first ``gap_after`` gaps in the
        nucleotide sequence and afterwards after each third gap,
        default is 2
    :param int tt: the number of the translation table, default is 1
    """
    gc = gcode(tt)
    aas = []
    ngap = 0
    check_start = check_start if check_start is not None else not complete
    codon = ''
    for i, nt in enumerate(str(seq).replace('U', 'T')):
        if nt == gap:
            ngap += 1
        else:
            codon = codon + nt
        if gap and gap_after is not None and ngap == gap_after:
            aas.append(gap)
            ngap -= 3
        if len(codon) == 3:
            if check_start:
                check_start = False
                if codon not in gc.starts and codon not in gc.astarts:
                    msg = (f'Codon {codon} is not a start codon {gc.starts} '
                           f'in genetic code #{gc.id}')
                    raise ValueError(msg)
                if codon not in gc.starts:
                    msg = f'Codon {codon} possibly is not a start codon.'
                    warnings.warn(msg)
            try:
                aa = gc.tt[codon]
            except (KeyError):
                aa = 'X'
            if codon in gc.astops:
                aa = astop
                if warn:
                    warnings.warn(f'Codon {codon} might be a stop codon.')
            if codon in gc.stops:
                if check_stop and i < len(nt) - 1:
                    msg = 'First stop codon is not at the end of the sequence.'
                    raise ValueError(msg)
                if not complete:
                    break
            aas.append(aa)
            codon = ''
    if check_stop and aa != 'X':
        msg = f'Last codon is not a stop codon {gc.stops} in genetic code #{gc.id}'
        raise ValueError(msg)
    return ''.join(aas)
