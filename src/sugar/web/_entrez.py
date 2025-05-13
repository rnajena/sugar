# (C) 2023, Tom Eulenfeld, MIT license
"""
Entrez client class

.. warning::
   This module is still experimental.
"""
import io
import os
from collections import deque
from time import perf_counter, sleep

from sugar import read, BioBasket


class Entrez():
    """
    Entrez client

    :param path: The path for persistence of downloaded files, default: no persistence,
        alternatively the path can be set with the ``ENTREZ_PAH`` environment variable.
    :param api_key: Optionally, you can use an API key that allows you to make
        more requests than without an API key,
        alternatively, the API key can be set with the ``ENTREZ_API_KEY`` environment variable.

    Without an API key you can make 3 requests per second,
    with an API key you can make 10 requests per second.
    The client will make sure that you do not use up this quota.
    By setting the path (environment) variable,
    repeated requests for the same id will not count against the quota.

    .. rubric:: Example:

    >>> from sugar.web import Entrez
    >>> client = Entrez()
    >>> seq = client.get_seq('AF522874')  # fetch multiple seqs with client.get_basket()

    """
    _url_fetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    _mail = 'tom.eulenfeld@uni-jena.de'
    _tool = 'sugar'
    _requests = 3
    _requests_api_key = 10
    _sleep = 1

    def __init__(self, path=os.getenv('ENTREZ_PATH'), api_key=os.getenv('ENTREZ_API_KEY')):
        self.path = path
        self.api_key = api_key
        self._request_times = deque()

    def _wait_before_request(self):
        """
        Ensure that not more requests than specified in the below URL are performed

        https://www.ncbi.nlm.nih.gov/books/NBK25497/
        """
        requests = self._requests_api_key if self.api_key else self._requests
        if len(self._request_times) >= requests:
            prev_time = self._request_times.popleft()
            time_elapsed = perf_counter() - prev_time
            if time_elapsed < self._sleep:
                sleep(self._sleep - time_elapsed)
        self._request_times.append(perf_counter())

    def fetch_seq(self, seqid, *, rettype='gb', db='nuccore', ext=None,
                  retmode='text', overwrite=False, path=None):
        r"""
        Fetch a sequence using the client

        :param seqid: Id of the sequence to be fetched

        :param str path: An alternative path for persistence,
            which might be different from the initialized path.
        :param bool overwrite: If True, redownload the sequence,
            even if it already exists in the path
        :param str ext: The file extension, defaults to rettype parameter.
        :param \*\*kw: Other kwargs are used to construct the request url,
            values other than the defaults are untested.

        :return: The filename of the downloaded content.
            If path is not set, the content is returned as ``StringIO`` instance.
        """
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
            self._wait_before_request()
            r = requests.get(self._url_fetch, params=params)
            r.raise_for_status()
            if fname is None:
                return io.StringIO(r.text)
            else:
                with open(fname, 'w') as f:
                    f.write(r.text)
        return fname

    def fetch_basket(self, seqids, **kw):
        r"""
        Fetch multiple sequences using the client

        :param seqids: A list of ids to fetch
        :param \*\*kw: All other kwargs are passed to `fetch_seq()`
        :returns: List of filenames or ``StringIO`` objects
        """
        return [self.fetch_seq(seqid, **kw) for seqid in seqids]

    def get_seq(self, seqid, *, read_kw=None, **kw):
        r"""
        Fetch a sequence and return it

        :param seqid: Id of the sequence to be fetched
        :param dict read_kw: Dictionary of read options passed to `.read()`
        :param \*\*kw: All other kwargs are passed to `fetch_seq()`
        :return: Fetched `.BioSeq` object
        """
        if read_kw is None:
            read_kw = dict()
        fname = self.fetch_seq(seqid, **kw)
        return read(fname, **read_kw)[0]

    def get_basket(self, seqids, *, read_kw=None, **kw):
        r"""
        Fetch multiple sequences and return them

        :param seqids: A list of ids to fetch
        :param dict read_kw: Dictionary of reading options passed to `.read()`
        :param \*\*kw: All other kwargs are passed to `fetch_basket()`
        :return: Fetched `.BioBasket` object
        """
        if read_kw is None:
            read_kw = dict()
        fnames = self.fetch_basket(seqids, **kw)
        return BioBasket([seq for fn in fnames for seq in read(fn, **read_kw)])
