"""Publication contracts using an in-memory VOSpace client."""

import csv
import hashlib
import io
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'workflow/scripts'))
from index_vos import reindex  # noqa: E402
from publish_vos import upload_once  # noqa: E402


class NotFound(Exception):
    pass


class FakeClient:
    def __init__(self):
        self.files = {}
        self.dirs = {'vos:root'}
        self.bad_md5 = False

    def get_node(self, uri, force=False):
        if uri in self.dirs:
            return SimpleNamespace(props={})
        if uri not in self.files:
            raise NotFound(uri)
        data = self.files[uri]
        return SimpleNamespace(props={'length': len(data), 'MD5': 'bad' if self.bad_md5 else hashlib.md5(data).hexdigest()})

    def mkdir(self, uri):
        self.dirs.add(uri)

    def copy(self, local, uri, send_md5=False):
        self.files[uri] = Path(local).read_bytes()

    def move(self, src, dst):
        self.files[dst] = self.files.pop(src)

    def listdir(self, uri, force=False):
        prefix = uri + '/'
        return sorted({path[len(prefix):].split('/')[0] for path in self.dirs | self.files.keys()
                       if path.startswith(prefix)})


def test_write_once_rejects_second_upload(tmp_path):
    client = FakeClient()
    local = tmp_path / 'pixel'
    local.write_bytes(b'first')
    uri = 'vos:root/catalogue/dataset/Norder=5/Dir=0/Npix=1.parquet'
    upload_once(client, local, uri)
    local.write_bytes(b'second')
    with pytest.raises(FileExistsError):
        upload_once(client, local, uri)
    assert client.files[uri] == b'first'


def test_md5_mismatch_never_exposes_final(tmp_path):
    client = FakeClient()
    client.bad_md5 = True
    local = tmp_path / 'pixel'
    local.write_bytes(b'contents')
    with pytest.raises(ValueError, match='MD5'):
        upload_once(client, local, 'vos:root/Npix=1.parquet')
    assert 'vos:root/Npix=1.parquet' not in client.files


def test_partition_list_from_listing_excludes_partial_and_local(tmp_path):
    client = FakeClient()
    root = 'vos:root/catalogue'
    directory = root + '/dataset/Norder=5/Dir=0'
    client.dirs.update({root, root + '/dataset', root + '/dataset/Norder=5', directory})
    client.files[directory + '/Npix=1.parquet'] = b'complete'
    client.files[directory + '/Npix=2.parquet.abc.part'] = b'partial'
    (tmp_path / 'Npix=3.parquet').write_bytes(b'not uploaded')
    assert reindex(client, root) == [1]
    rows = list(csv.DictReader(io.StringIO(client.files[root + '/partition_info.csv'].decode())))
    assert rows == [{'Norder': '5', 'Npix': '1'}]


def test_owner_rejects_other_and_duplicate_claims(tmp_path):
    from publish_pixel import require_owner
    claims = tmp_path / 'claims.yaml'
    claims.write_text('order: 2\ncells:\n  1: nibi\n')
    require_owner(64, 'nibi', claims)
    with pytest.raises(ValueError):
        require_owner(64, 'fir', claims)
    claims.write_text('order: 2\ncells:\n  1: nibi\n  1: fir\n')
    with pytest.raises(ValueError, match='Duplicate'):
        require_owner(64, 'nibi', claims)
