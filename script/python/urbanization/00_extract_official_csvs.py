#!/usr/bin/env python3
                                                                              

                                                                           
                                                                         
                                                                           
   

from __future__ import annotations

import binascii
import concurrent.futures
import struct
import urllib.request
import zlib
from pathlib import Path


ARCHIVES = {
    "general": "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_UCDB_GLOBE_R2024A/GHS_UCDB_THEME_GLOBE_R2024A/GHS_UCDB_THEME_GENERAL_CHARACTERISTICS_GLOBE_R2024A/V1-2/GHS_UCDB_THEME_GENERAL_CHARACTERISTICS_GLOBE_R2024A_V1_2.zip",
    "ghsl": "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_UCDB_GLOBE_R2024A/GHS_UCDB_THEME_GLOBE_R2024A/GHS_UCDB_THEME_GHSL_GLOBE_R2024A/V1-2/GHS_UCDB_THEME_GHSL_GLOBE_R2024A_V1_2.zip",
    "greenness": "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_UCDB_GLOBE_R2024A/GHS_UCDB_THEME_GLOBE_R2024A/GHS_UCDB_THEME_GREENNESS_GLOBE_R2024A/V1-2/GHS_UCDB_THEME_GREENNESS_GLOBE_R2024A_V1_2.zip",
    "infrastructure": "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_UCDB_GLOBE_R2024A/GHS_UCDB_THEME_GLOBE_R2024A/GHS_UCDB_THEME_INFRASTRUCTURES_GLOBE_R2024A/V1-2/GHS_UCDB_THEME_INFRASTRUCTURES_GLOBE_R2024A_V1_2.zip",
}


def request(url: str, start: int | None = None, end: int | None = None, method: str = "GET"):
    headers = {"User-Agent": "urban-bird-virome/1.0"}
    if start is not None:
        headers["Range"] = f"bytes={start}-{'' if end is None else end}"
    req = urllib.request.Request(url, headers=headers, method=method)
    return urllib.request.urlopen(req, timeout=120)


def content_length(url: str) -> int:
    with request(url, method="HEAD") as response:
        if response.headers.get("Accept-Ranges", "").lower() != "bytes":
            raise RuntimeError(f"Server does not advertise byte ranges: {url}")
        return int(response.headers["Content-Length"])


def get_range(url: str, start: int, end: int) -> bytes:
    with request(url, start, end) as response:
        data = response.read()
    expected = end - start + 1
    if len(data) != expected:
        raise RuntimeError(f"Range length mismatch: expected {expected}, received {len(data)}")
    return data


def get_range_parallel(url: str, start: int, end: int, block_size: int = 1024 * 1024) -> bytes:
    ranges = []
    position = start
    while position <= end:
        block_end = min(end, position + block_size - 1)
        ranges.append((position, block_end))
        position = block_end + 1
    if len(ranges) == 1:
        return get_range(url, start, end)
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(8, len(ranges))) as executor:
        futures = [executor.submit(get_range, url, block_start, block_end) for block_start, block_end in ranges]
        blocks = [future.result() for future in futures]
    return b"".join(blocks)


def central_directory(url: str, size: int) -> bytes:
    tail_start = max(0, size - 131072)
    tail = get_range(url, tail_start, size - 1)
    signature = b"PK\x05\x06"
    position = tail.rfind(signature)
    if position < 0:
        raise RuntimeError("ZIP end-of-central-directory record not found")
    eocd = tail[position : position + 22]
    _, _, _, _, _, directory_size, directory_offset, _ = struct.unpack("<4s4H2IH", eocd)
    if directory_offset >= tail_start and directory_offset + directory_size <= size:
        local = directory_offset - tail_start
        return tail[local : local + directory_size]
    return get_range(url, directory_offset, directory_offset + directory_size - 1)


def csv_entry(directory: bytes) -> dict[str, int | str]:
    offset = 0
    entries = []
    while offset < len(directory):
        header = directory[offset : offset + 46]
        if len(header) != 46 or header[:4] != b"PK\x01\x02":
            raise RuntimeError(f"Invalid central-directory header at byte {offset}")
        values = struct.unpack("<4s6H3I5H2I", header)
        compression = values[4]
        crc32 = values[7]
        compressed_size = values[8]
        uncompressed_size = values[9]
        filename_length = values[10]
        extra_length = values[11]
        comment_length = values[12]
        local_offset = values[16]
        start = offset + 46
        filename_bytes = directory[start : start + filename_length]
        filename = filename_bytes.decode("utf-8")
        entries.append(
            {
                "filename": filename,
                "compression": compression,
                "crc32": crc32,
                "compressed_size": compressed_size,
                "uncompressed_size": uncompressed_size,
                "local_offset": local_offset,
            }
        )
        offset = start + filename_length + extra_length + comment_length
    matches = [entry for entry in entries if str(entry["filename"]).lower().endswith(".csv")]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one CSV member, found {[entry['filename'] for entry in matches]}")
    return matches[0]


def extract_csv(url: str, output_directory: Path) -> Path:
    size = content_length(url)
    entry = csv_entry(central_directory(url, size))
    local_offset = int(entry["local_offset"])
    local_header = get_range(url, local_offset, local_offset + 29)
    values = struct.unpack("<4s5H3I2H", local_header)
    if values[0] != b"PK\x03\x04":
        raise RuntimeError("Invalid local ZIP header")
    filename_length = values[9]
    extra_length = values[10]
    data_start = local_offset + 30 + filename_length + extra_length
    compressed_size = int(entry["compressed_size"])
    compressed = get_range_parallel(url, data_start, data_start + compressed_size - 1)
    compression = int(entry["compression"])
    if compression == 0:
        raw = compressed
    elif compression == 8:
        raw = zlib.decompress(compressed, -zlib.MAX_WBITS)
    else:
        raise RuntimeError(f"Unsupported ZIP compression method: {compression}")
    if len(raw) != int(entry["uncompressed_size"]):
        raise RuntimeError("Uncompressed CSV size does not match ZIP directory")
    crc = binascii.crc32(raw) & 0xFFFFFFFF
    if crc != int(entry["crc32"]):
        raise RuntimeError("CSV CRC check failed")
    output = output_directory / str(entry["filename"])
    output.write_bytes(raw)
    return output


def main() -> None:
    output_directory = Path(__file__).resolve().parent / "official_source_csv"
    output_directory.mkdir(parents=True, exist_ok=True)
    for name, url in ARCHIVES.items():
        extract_csv(url, output_directory)


if __name__ == "__main__":
    main()
