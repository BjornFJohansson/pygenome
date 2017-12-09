import requests


def do_thing():
    return requests.get('https://httpbin.org/get').json()