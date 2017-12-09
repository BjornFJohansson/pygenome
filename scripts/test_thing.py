from . import thing

import pytest
import requests_mock as rm_module


@pytest.fixture
def requests_mock(request):
    m = rm_module.Mocker()
    m.start()
    request.addfinalizer(m.stop)
    return m


def test_a_thing(requests_mock):
    requests_mock.get('https://httpbin.org/get', json={'last-modified': 'thing'})
    json = thing.do_thing()
    assert json == {'fake': 'thing'}
    print(123)