import pytest


def pytest_addoption(parser):
    parser.addoption('--web', action='store_true',
                 default=False, help='enable web tests')

def pytest_collection_modifyitems(config, items):
    if not config.getoption('--web'):
        skip_web = pytest.mark.skip(reason='webtest: need --web option to run')
        for item in items:
            if 'webtest' in item.keywords:
                item.add_marker(skip_web)
    # explicitely add filter warnings to markers so that they have a higher
    # priority than command line options, e.g. -W error
    for item in items:
        for fwarn in config.getini('filterwarnings'):
            item.add_marker(pytest.mark.filterwarnings(fwarn))
