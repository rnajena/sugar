# (C) 2024, Tom Eulenfeld, MIT license
"""
Home of the ``sugar`` command line script

Run ``sugar -h`` to see all available subcommands.
``sugar test`` runs the test suite.
"""

import argparse
import sys
import contextlib
from pathlib import Path
import os


@contextlib.contextmanager
def _changedir(path):
    origin = Path().resolve()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)


def _start_ipy(seqs):
    from IPython import start_ipython
    print('Contents loaded into seqs variable.')
    start_ipython(argv=[], user_ns={'seqs': seqs}, display_banner=False)
    print('Bye')


def cat(fname, fmt=None, out=None, fmtout=None, merge=False):
    """Concatenate sequence files, optionally merge and convert to different format (sugar cat, sugar merge)"""
    from sugar import BioBasket, read
    seqs = BioBasket()
    for fn in fname:
        seqs += read(fn, fmt)
    if merge:
        seqs.merge()
    if out is None:
        try:
            print(seqs.tofmtstr(fmtout or fmt or seqs[0].meta._fmt))
        except BrokenPipeError:
            pass
    else:
        seqs.write(out, fmt=fmtout)


def catf(fname, fmt=None, out=None, fmtout=None):
    """Concatenate fts files and convert to different format (sugar catf)"""
    from sugar import FeatureList, read_fts
    fts = FeatureList()
    for fn in fname:
        fts += read_fts(fn, fmt)
    if out is None:
        try:
            print(fts.tofmtstr(fmtout or fmt or fts[0].meta._fmt))
        except BrokenPipeError:
            pass
    else:
        fts.write(out, fmt=fmtout)


def translate(fname, fmt, out=None, fmtout=None, cds=False, **kw):
    from sugar import read

    try:
        seqs = read(fname, fmt)
    except Exception as ex1:
        try:
            from sugar.core.cane import translate as trans
            for line in fname.splitlines():
                print(trans(line, **kw))
            return
        except Exception as ex2:
            raise ExceptionGroup(
                'Expect sequence file or nucleotide string as input',
                [ex1, ex2])
    if cds:
        seqs = seqs['cds']
    seqs.translate(**kw)
    if out is None:
        try:
            print(seqs.tofmtstr(fmtout or fmt or seqs[0].meta._fmt))
        except BrokenPipeError:
            pass
    else:
        seqs.write(out, fmt=fmtout)


def run(command, pytest_args=None, pdb=False, fname=None, fmt=None, **kw):
    """Dispatch command to function"""
    if pdb:
        import traceback, pdb
        def info(type, value, tb):
            traceback.print_exception(type, value, tb)
            print()
            pdb.pm()
        sys.excepthook = info
    if command == 'print':
        from sugar import read
        try:
            print(read(fname, fmt).tostr(**kw))
        except BrokenPipeError:
            pass
    elif command == 'printf':
        from sugar import read_fts
        try:
            print(read_fts(fname, fmt).tostr(**kw))
        except BrokenPipeError:
            pass
    elif command == 'load':
        from sugar import read
        seqs = read(fname, fmt, **kw)
        _start_ipy(seqs)
    elif command == 'loadf':
        from sugar import read_fts
        fts = read_fts(fname, fmt, **kw)
        _start_ipy(fts)
    elif command in ('cat', 'merge'):
        cat(fname, fmt=fmt, merge=(command=='merge'), **kw)
    elif command == 'catf':
        catf(fname, fmt=fmt, **kw)
    elif command == 'index':
        from sugar.index.fastaindex import _fastaindex_cmd
        _fastaindex_cmd(**kw)
    elif command == 'translate':
        translate(fname, fmt=fmt, **kw)
    elif command == 'test':
        try:
            import pytest
        except ImportError:
            msg = ("\nsugar's test suite uses pytest. "
                    "Please install pytest before running the tests.")
            sys.exit(msg)
        path = Path(__file__).parent / 'tests'
        print(f'Run pytest in directory {path}')
        with _changedir(path):
            status = pytest.main(pytest_args)
        sys.exit(status)
    else:
        raise ValueError(f'Unknown command: {command}')


def _str2tuple(t):
    return tuple(tt.strip() for tt in t.split(','))


def cli(cmd_args=None):
    """Main entry point from the command line"""
    from sugar import __version__
    msg = ('Sugar for your RNA')
    epilog = 'To get help on a subcommand run: sugar command -h'
    parser = argparse.ArgumentParser(description=msg, epilog=epilog)
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('--pdb', action='store_true', help='Start the debugger upon exception')

    sub = parser.add_subparsers(title='commands', dest='command')
    sub.required = True
    p_cat = sub.add_parser('cat', help='concatenate sequence files or convert to different format')
    p_catf = sub.add_parser('catf', help='concatenate fts files or convert to different format')
    p_index = sub.add_parser('index', help='index FASTA files, query the database')
    p_load = sub.add_parser('load', help='load seq file into IPython session')
    p_loadf = sub.add_parser('loadf', help='load fts file into IPython session')
    p_merge = sub.add_parser('merge', help='merge sequences with the same id')
    p_print = sub.add_parser('print', help='print contents of seq file')
    p_printf = sub.add_parser('printf', help='print contents of fts file')
    p_trans = sub.add_parser('translate', help='translate nucleotide sequence')
    msg = 'run sugar test suite'
    msg2 = ('The test suite uses pytest. You can call pytest directly or use '
            'most of pytest cli arguments in the sugar test call. '
            'See pytest -h. Use the --web option to additionally run web tests.')
    p_test = sub.add_parser('test', help=msg, description=msg2)

    for p in (p_print, p_printf, p_load, p_loadf):
        p.add_argument('fname', help='filenames')
        p.add_argument('-f', '--fmt', help='format supported by Sugar (default: auto-detect)')
    for p in (p_print, p_load):
        p.add_argument('-e', '--exclude', help='exclude contents', type=_str2tuple, default=argparse.SUPPRESS)
    p_print.add_argument('--raw', help='just print the letters, one sequence on each line', action='store_true')
    p_printf.add_argument('--raw', help='just print a table with type, start index, stop index, name', action='store_true')
    for p in (p_print, p_printf):
        p.add_argument('--h', help='max height, 0 for no restriction, default 19', default=19, type=int)
        p.add_argument('--w', help='max width, 0 for no restriction, default 80', default=80, type=int)
    p_print.add_argument('--wid', help='max id width, 0 for no restriction, default 19', default=19, type=int)
    # p_print.add_argument('--showgc', help='show GC content', default=True, action=argparse.BooleanOptionalAction)
    p_print.add_argument('--no-showgc', help='do not show GC content', action='store_false', dest='showgc')
    for p in (p_cat, p_catf, p_merge):
        p.add_argument('fname', help='filename(s) in', nargs='+')
    p_trans.add_argument('fname', help='filename in')
    for p in (p_cat, p_catf, p_merge, p_trans):
        p.add_argument('-o', '--out', help='filename out')
        p.add_argument('-f', '--fmt', help='format in')
        p.add_argument('-fo', '--fmtout', help='format out')

    p_trans.add_argument('-tt', '--translation-table', help='number of translation table, default 1', default=1, type=int, dest='tt')
    p_trans.add_argument('-c', '--complete', help='whether to ignore stop codons', action='store_true')
    p_trans.add_argument('--cds', help='cut out CDS feature before translation', action='store_true')

    sub_index = p_index.add_subparsers(title='commands', dest='idxcommand')
    p_idxinfo = sub_index.add_parser('info', help='show info about index file')
    p_idxcreate = sub_index.add_parser('create', help='create new index file')
    p_idxcreate.add_argument('-m', '--mode', choices=('binary', 'db'), default='binary')
    p_idxcreate.add_argument('-p', '--path', default='{dbpath}')
    p_idxadd = sub_index.add_parser('add', help='add fasta files')
    p_idxadd.add_argument('fnames', nargs='+')
    p_idxfetch = sub_index.add_parser('fetch', help='fetch fasta')
    p_idxprint = sub_index.add_parser('print', help='print seqs')
    p_idxload = sub_index.add_parser('load', help='load seqs')
    for p in (p_idxfetch, p_idxprint, p_idxload):
        p.add_argument('seqids')
    p_idxcreate.add_argument('dbname', help='database file')
    for p in (p_idxinfo, p_idxadd, p_idxfetch, p_idxprint, p_idxload):
        p.add_argument('-d', '--dbname', help='database file, by default last used')
    p_idxfetch.add_argument('-o', '--out', default='-')

    # Get command line arguments and start run function
    args, pytest_args = parser.parse_known_args(cmd_args)
    if args.command != 'test':
        parser.parse_args(cmd_args)  # just call again to properly raise the error
    run(pytest_args=pytest_args, **vars(args))


if __name__ == '__main__':
    cli()
