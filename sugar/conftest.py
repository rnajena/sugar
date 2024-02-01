import pytest


def pytest_addoption(parser):
    parser.addoption('--web', action='store_true',
                 default=False, help='enable web tests')

def pytest_configure(config):
    config.addinivalue_line('markers', 'webtest: mark test as a web test')

def pytest_collection_modifyitems(config, items):
    if not config.getoption('--web'):
        skip_web = pytest.mark.skip(reason='webtest: need --web option to run')
        for item in items:
            if 'webtest' in item.keywords:
                item.add_marker(skip_web)
