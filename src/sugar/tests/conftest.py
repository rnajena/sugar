import pytest
from copy import deepcopy


def pytest_addoption(parser):
    parser.addoption('--web', action='store_true',
                 default=False, help='enable web tests')

def pytest_collection_modifyitems(config, items):
    if not config.getoption('--web'):
        skip_web = pytest.mark.skip(reason='webtest: require --web option')
        for item in items:
            if 'webtest' in item.keywords:
                item.add_marker(skip_web)
    # explicitely add filter warnings to markers so that they have a higher
    # priority than command line options, e.g. -W error
    for item in items:
        for fwarn in config.getini('filterwarnings'):
            item.add_marker(pytest.mark.filterwarnings(fwarn))


@pytest.fixture(scope='session', autouse=True)
def fix_entrypoints_in_pytest_session():
    from importlib.metadata import EntryPoint, EntryPoints
    from sugar._io import util
    ep = EntryPoint(name='testbin', value='sugar.tests.test_io_binary_plugin', group='sugar.io')
    ep2 = EntryPoint(name='testbin', value='sugar.tests.test_io_binary_plugin', group='sugar.io.fts')
    ep3 = EntryPoint(name='fancy', value='sugar.tests.test_io_template_plugin', group='sugar.io')
    original_EPS = deepcopy(util.EPS)
    util.EPS['seqs'] = EntryPoints(util.EPS['seqs'] + (ep, ep3))
    util.EPS['fts'] = EntryPoints(util.EPS['fts'] + (ep2,))
    util.FMTS['seqs'].extend(['testbin', 'fancy'])
    util.FMTS_ALL['seqs'].extend(['testbin', 'fancy'])
    util.FMTS['fts'].append('testbin')
    util.FMTS_ALL['fts'].append('testbin')
    yield
    util.EPS = original_EPS
    util.FMTS['seqs'].remove('testbin')
    util.FMTS_ALL['seqs'].remove('testbin')
    util.FMTS['fts'].remove('testbin')
    util.FMTS_ALL['fts'].remove('testbin')
