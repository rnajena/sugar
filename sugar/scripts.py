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


def convert(fname, fmt=None, out=None, fmtout=None,
            tool=None, toolout=None):
    """Convert files to different format (sugar convert)"""
    from sugar import read
    seqs = read(fname, fmt, tool=tool)
    if out is None and fmtout is None:
        raise ValueError('fmt need to be specified for out=None')
    elif out is None:
        try:
            print(seqs.tofmtstr(fmtout))
        except BrokenPipeError:
            pass
    else:
        seqs.write(out, fmt=fmtout, tool=toolout)


def run(command, pytest_args=None, pdb=False, fname=None, fmt=None, tool=None, **kw):
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
            print(read(fname, fmt, tool=tool).tostr(**kw))
        except BrokenPipeError:
            pass
    elif command == 'load':
        from sugar import read
        seqs = read(fname, fmt, tool=tool, **kw)
        _start_ipy(seqs)
    elif command == 'convert':
        convert(fname, fmt=fmt, tool=tool, **kw)
    elif command == 'index':
        from sugar._io.fastaindex import _fastaindex_cmd
        _fastaindex_cmd(**kw)
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


def _str2tuple(t):
    return tuple(tt.strip() for tt in t.split(','))


def cli(cmd_args=None):
    """Main entry point from the command line"""
    from sugar import __version__
    msg = ('Sugar for your RNA')
    epilog = 'To get help on a subcommand run: sugar command -h'
    parser = argparse.ArgumentParser(description=msg, epilog=epilog)
    version = '%(prog)s ' + __version__
    parser.add_argument('--version', action='version', version=version)
    msg = 'Start the debugger upon exception'
    parser.add_argument('--pdb', action='store_true', help=msg)

    sub = parser.add_subparsers(title='commands', dest='command')
    sub.required = True
    msg = 'print contents of seq file'
    p_print = sub.add_parser('print', help=msg)
    msg = 'convert between different file formats'
    p_convert = sub.add_parser('convert', help=msg)
    msg = 'load objects into IPython session'
    p_load = sub.add_parser('load', help=msg)
    p_index = sub.add_parser('index', help='index FASTA files, query the database')
    msg = 'run sugar test suite'
    msg2 = ('The test suite uses pytest. You can call pytest directly or use '
            'most of pytest cli arguments in the sugar test call. '
            'See pytest -h. Use the --web option to additionally run web tests.')
    p_test = sub.add_parser('test', help=msg, description=msg2)

    for p in (p_print, p_load):
        p.add_argument('fname', help='filenames')
        msg = 'format supported by Sugar (default: auto-detect)'
        p.add_argument('-f', '--fmt', help=msg)
        p.add_argument('-t', '--tool', help='tool for reading')
        p.add_argument('-e', '--exclude', help='exclude contents', type=_str2tuple, default=argparse.SUPPRESS)
    p_print.add_argument('--h', help='max height, 0 for no restriction, default 19', default=19, type=int)
    p_print.add_argument('--w', help='max width, 0 for no restriction, default 80', default=80, type=int)
    p_print.add_argument('--wid', help='max id width, 0 for no restriction, default 19', default=19, type=int)
    # p_print.add_argument('--showgc', help='show GC content', default=True, action=argparse.BooleanOptionalAction)
    p_print.add_argument('--no-showgc', help='do not show GC content', action='store_false', dest='showgc')
    p_convert.add_argument('fname', help='filename in')
    p_convert.add_argument('-o', '--out', help='filename out')
    p_convert.add_argument('-f', '--fmt', help='format in')
    p_convert.add_argument('-fo', '--fmtout', help='format out')
    p_convert.add_argument('-t', '--tool', help='tool for reading')
    p_convert.add_argument('-to', '--toolout', help='tool for writing')


    sub_index = p_index.add_subparsers(title='commands', dest='idxcommand')
    p_idxinfo = sub_index.add_parser('info')
    p_idxcreate = sub_index.add_parser('create', help='create stuff')
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
