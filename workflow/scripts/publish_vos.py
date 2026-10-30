#!/usr/bin/env python3
"""Write immutable publication files to VOSpace using an injected vos client."""

import hashlib
import time
from pathlib import Path
from uuid import uuid4


def exists(client, uri):
    """Return whether a node exists; propagate errors other than not-found."""
    try:
        client.get_node(uri, force=True)
        return True
    except Exception as exc:
        if '404' in str(exc) or 'NotFound' in type(exc).__name__:
            return False
        raise


def ensure_dir(client, uri):
    """Create a missing container and its ancestors."""
    if exists(client, uri):
        return
    parent = uri.rsplit('/', 1)[0]
    if parent != uri:
        ensure_dir(client, parent)
    try:
        client.mkdir(uri)
    except FileExistsError:
        # Another pixel publisher may have created the shared container.
        if not exists(client, uri):
            raise


def verify(client, local, uri):
    """Require matching size and MD5, allowing bounded CADC property lag.

    Re-read the node up to 13 times, with at most 60 seconds of retry waits.
    A persistent mismatch still fails closed; no final name is exposed before
    the temporary node verifies. Hash the local file only once per call.
    """
    digest = hashlib.md5()  # VOSpace exposes MD5 for transfer verification.
    with Path(local).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    length, md5 = Path(local).stat().st_size, digest.hexdigest()
    for attempt in range(13):
        props = client.get_node(uri, force=True).props
        remote_md5 = str(props.get('MD5') or '').lower().removeprefix('md5:')
        if int(props.get('length') or 0) == length and remote_md5 == md5:
            return
        if attempt < 12:
            time.sleep(5)
    raise ValueError(f'VOSpace size or MD5 mismatch: {uri}')


def upload_once(client, local, uri):
    """Upload and verify a temporary node before exposing its final name.

    @sc [label:publication] write-once
    Pixel and map destinations are immutable: an existing node is refused,
    and the uploaded bytes must match in size and MD5 before and after move.
    tests/unit/test_publish_vos.py checks refusal and checksum mismatch.
    """
    if exists(client, uri):
        raise FileExistsError(uri)
    ensure_dir(client, uri.rsplit('/', 1)[0])
    temporary = uri + '.' + uuid4().hex + '.part'
    client.copy(str(local), temporary, send_md5=True)
    verify(client, local, temporary)
    if exists(client, uri):
        raise FileExistsError(uri)
    client.move(temporary, uri)
    verify(client, local, uri)
