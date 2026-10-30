#!/usr/bin/env python3
"""Rebuild a HATS partition list from VOSpace's published leaf listing."""

import csv
import re
import tempfile
from pathlib import Path

from publish_vos import exists


def partitions(client, root):
    """List complete order-5 HATS leaves, ignoring temporary names.

    @sc [label:publication] partition-list-from-listing
    The partition list is derived solely from remote final leaf names, never
    local builds or partial uploads. tests/unit/test_publish_vos.py exercises
    a listing with temporary and locally present but unuploaded files.
    """
    dataset = root.rstrip('/') + '/dataset/Norder=5'
    if not exists(client, dataset):
        return []
    pixels = set()
    for directory in client.listdir(dataset, force=True):
        directory = directory.rstrip('/').rsplit('/', 1)[-1]
        match = re.fullmatch(r'Dir=(\d+)', directory)
        if not match:
            continue
        for name in client.listdir(dataset + '/' + directory, force=True):
            name = name.rsplit('/', 1)[-1]
            leaf = re.fullmatch(r'Npix=(\d+)\.parquet', name)
            if leaf:
                pixel = int(leaf.group(1))
                if 0 <= pixel < 12 * 4 ** 5 and pixel // 10000 * 10000 == int(match.group(1)):
                    pixels.add(pixel)
    return sorted(pixels)


def reindex(client, root):
    """Overwrite only partition_info.csv with the remote listing."""
    pixels = partitions(client, root)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / 'partition_info.csv'
        with path.open('w', newline='') as stream:
            writer = csv.writer(stream, lineterminator='\n')
            writer.writerow(['Norder', 'Npix'])
            writer.writerows((5, pixel) for pixel in pixels)
        client.copy(str(path), root.rstrip('/') + '/partition_info.csv', send_md5=True)
    return pixels
