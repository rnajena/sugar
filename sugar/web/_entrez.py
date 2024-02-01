# (C) 2023, Tom Eulenfeld, MIT license
import io
import os
from collections import deque
from time import perf_counter, sleep

from sugar import read, BioBasket


class Entrez():
    _url_fetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    _mail = 'tom.eulenfeld@uni-jena.de'
    _tool = 'sugar'
    _requests = 3
    _requests_api_key = 10
    _seconds = 1

    def __init__(self, path=os.getenv('ENTREZ_PATH'), api_key=os.getenv('ENTREZ_API_KEY')):
        self.path = path
        self.api_key = api_key
        self._request_times = deque()

    def wait_before_request(self):
        requests = self._requests_api_key if self.api_key else self._requests
        if len(self._request_times) >= requests:
            prev_time = self._request_times.popleft()
            time_elapsed = perf_counter() - prev_time
            if time_elapsed < self._seconds:
                sleep(self._seconds - time_elapsed)
        self._request_times.append(perf_counter())

    def fetch_seq(self, seqid, *, rettype='gb', db='nuccore', ext=None,
                  retmode='text', overwrite=False, path=None):
        import requests

        path = path or self.path
        if ext is None:
            ext = rettype
        if path and not os.path.isdir(str(path)):
            os.makedirs(str(path))
        if path is not None:
            fname = os.path.join(path, seqid + '.' + ext)
        else:
            fname = None
        if (fname is None or not os.path.isfile(fname) or
                os.path.getsize(fname) == 0 or overwrite):
            params = dict(db=db, id=seqid, rettype=rettype, retmode=retmode,
                          mail=self._mail, tool=self._tool)
            if self.api_key:
                params['api_key'] = self.api_key
            self.wait_before_request()
            r = requests.get(self._url_fetch, params=params)
            r.raise_for_status()
            if fname is None:
                return io.StringIO(r.text)
            else:
                with open(fname, 'w') as f:
                    f.write(r.text)
        return fname

    def fetch_basket(self, seqids, **kw):
        return [self.fetch_seq(seqid, **kw) for seqid in seqids]

    def get_seq(self, seqid, *, read_kw=None, **kw):
        if read_kw is None:
            read_kw = dict()
        fname = self.fetch_seq(seqid, **kw)
        return read(fname, **read_kw)[0]

    def get_basket(self, seqids, *, read_kw=None, **kw):
        if read_kw is None:
            read_kw = dict()
        fnames = self.fetch_basket(seqids, **kw)
        return BioBasket([seq for fn in fnames for seq in read(fn, **read_kw)])
