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


@pytest.fixture(autouse=True)
def no_verification_wait(monkeypatch):
    """Keep bounded retry tests fast without weakening their assertions."""
    monkeypatch.setattr('publish_vos.time.sleep', lambda seconds: None)


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


def test_property_lag_before_and_after_move(tmp_path, monkeypatch):
    class DelayedClient(FakeClient):
        def __init__(self):
            super().__init__()
            self.reads = {}

        def get_node(self, uri, force=False):
            node = super().get_node(uri, force=force)
            if uri in self.files:
                self.reads[uri] = self.reads.get(uri, 0) + 1
                if self.reads[uri] <= 2:
                    return SimpleNamespace(props={'length': 0, 'MD5': ''})
            return node

    sleeps = []
    monkeypatch.setattr('publish_vos.time.sleep', sleeps.append)
    client = DelayedClient()
    local = tmp_path / 'pixel'
    local.write_bytes(b'contents')
    upload_once(client, local, 'vos:root/Npix=1.parquet')
    assert client.files == {'vos:root/Npix=1.parquet': b'contents'}
    assert sleeps == [5, 5, 5, 5]
    assert sorted(client.reads.values()) == [3, 3]


def test_property_lag_is_bounded(tmp_path, monkeypatch):
    sleeps = []
    monkeypatch.setattr('publish_vos.time.sleep', sleeps.append)
    client = FakeClient()
    client.bad_md5 = True
    local = tmp_path / 'pixel'
    local.write_bytes(b'contents')
    with pytest.raises(ValueError, match='MD5'):
        upload_once(client, local, 'vos:root/Npix=1.parquet')
    assert sleeps == [5] * 12
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


def test_publish_chains_build_upload_index_and_maps(tmp_path, monkeypatch):
    import publish_pixel

    client = FakeClient()
    claims = tmp_path / 'claims.yaml'
    claims.write_text('order: 2\ncells:\n  0: nibi\n')
    args = SimpleNamespace(pixel=1, cluster='nibi', claims=claims,
                           output=tmp_path / 'out', root='vos:root',
                           tile_list='tiles', exposure_list='exposures',
                           products_dir='products', run_dir='run', sha='abc', mask=['mask'],
                           galaxy_schema_from=None, star_schema_from=None)

    def build(pixel, kind, source, products, run, output, sha, mask, schema_from=None):
        path = output / 'dataset/Norder=5/Dir=0/Npix=1.parquet'
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(kind.encode())
        return path

    def metadata(root, sample, kind, name):
        root = Path(root)
        (root / 'dataset').mkdir()
        (root / 'properties').write_text(kind)
        (root / 'dataset/_common_metadata').write_text(kind)

    def maps(pixel, source, products, run, output, sha, mask):
        paths = {}
        for name in ('npoint', 'ntaken', 'nexp', 'nflagged'):
            path = output / name / 'Npix=1.hsp'
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(name)
            paths[name] = path
        return paths

    monkeypatch.setattr(publish_pixel, 'build_pixel', build)
    monkeypatch.setattr(publish_pixel, 'write_catalogue_metadata', metadata)
    monkeypatch.setattr(publish_pixel, 'build_maps', maps)
    publish_pixel.publish(client, args)
    for name in ('galaxies', 'stars'):
        assert b'5,1' in client.files[f'vos:root/{name}/partition_info.csv']
    assert client.files['vos:root/maps/nexp/Npix=1.hsp'] == b'nexp'


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
