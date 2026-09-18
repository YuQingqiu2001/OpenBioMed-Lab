"""Physical-coordinate slide readers for SVS and SR386 CZI.

All public coordinates are micrometres from a normalized top-left slide origin.
The native CZI directory may use negative starts; these are normalized once and
recorded in the manifest instead of leaking into downstream topology.
"""

from __future__ import annotations

import math
import sqlite3
import xml.etree.ElementTree as ET
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import cv2
import numpy as np
from PIL import Image


@dataclass(frozen=True)
class SlideGeometry:
    path: str
    reader: str
    native_mpp_x: float
    native_mpp_y: float
    native_origin_x_px: int
    native_origin_y_px: int
    width_px: int
    height_px: int

    @property
    def width_um(self) -> float:
        return self.width_px * self.native_mpp_x

    @property
    def height_um(self) -> float:
        return self.height_px * self.native_mpp_y

    def to_dict(self) -> Dict[str, Any]:
        result = dict(self.__dict__)
        result["width_um"] = self.width_um
        result["height_um"] = self.height_um
        return result


class PhysicalSlideReader:
    geometry: SlideGeometry

    def read_bbox_um(
        self, bbox_um: Tuple[float, float, float, float], target_mpp: float
    ) -> np.ndarray:
        raise NotImplementedError

    def coarse_tissue_presence(self, target_mpp: float = 8.0) -> np.ndarray:
        """Return a conservative low-resolution tissue/presence mask."""

        rgb = self.read_bbox_um(
            (0.0, 0.0, self.geometry.width_um, self.geometry.height_um), target_mpp
        )
        array = rgb.astype(np.float32)
        brightness = array.mean(axis=2)
        chroma = array.max(axis=2) - array.min(axis=2)
        return (brightness < 238.0) & (chroma > 6.0)

    def close(self) -> None:
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()


class OpenSlidePhysicalReader(PhysicalSlideReader):
    def __init__(self, path: Path, mpp_x: float | None = None, mpp_y: float | None = None):
        import openslide

        self.path = Path(path)
        self.slide = openslide.OpenSlide(str(path))
        props = self.slide.properties
        resolved_x = float(mpp_x or props.get(openslide.PROPERTY_NAME_MPP_X, "nan"))
        resolved_y = float(mpp_y or props.get(openslide.PROPERTY_NAME_MPP_Y, resolved_x))
        if not np.isfinite(resolved_x) or not np.isfinite(resolved_y):
            raise ValueError(f"Missing physical MPP for {path}")
        width, height = self.slide.dimensions
        self.geometry = SlideGeometry(
            path=str(path),
            reader="openslide",
            native_mpp_x=resolved_x,
            native_mpp_y=resolved_y,
            native_origin_x_px=0,
            native_origin_y_px=0,
            width_px=int(width),
            height_px=int(height),
        )

    def read_bbox_um(
        self, bbox_um: Tuple[float, float, float, float], target_mpp: float
    ) -> np.ndarray:
        x0_um, y0_um, x1_um, y1_um = bbox_um
        x0_um = max(0.0, x0_um)
        y0_um = max(0.0, y0_um)
        x1_um = min(self.geometry.width_um, x1_um)
        y1_um = min(self.geometry.height_um, y1_um)
        out_w = max(1, int(round((x1_um - x0_um) / target_mpp)))
        out_h = max(1, int(round((y1_um - y0_um) / target_mpp)))

        source_mpp = max(self.geometry.native_mpp_x, self.geometry.native_mpp_y)
        desired_downsample = max(1.0, target_mpp / source_mpp)
        level = self.slide.get_best_level_for_downsample(desired_downsample)
        level_downsample = float(self.slide.level_downsamples[level])
        level_w = max(
            1,
            int(
                math.ceil(
                    (x1_um - x0_um)
                    / (self.geometry.native_mpp_x * level_downsample)
                )
            ),
        )
        level_h = max(
            1,
            int(
                math.ceil(
                    (y1_um - y0_um)
                    / (self.geometry.native_mpp_y * level_downsample)
                )
            ),
        )
        x0_px = int(math.floor(x0_um / self.geometry.native_mpp_x))
        y0_px = int(math.floor(y0_um / self.geometry.native_mpp_y))
        rgba = self.slide.read_region((x0_px, y0_px), level, (level_w, level_h))
        rgb = np.asarray(rgba.convert("RGB"), dtype=np.uint8)
        if rgb.shape[:2] != (out_h, out_w):
            rgb = cv2.resize(rgb, (out_w, out_h), interpolation=cv2.INTER_AREA)
        return rgb

    def close(self) -> None:
        self.slide.close()


class IblPhysicalReader(PhysicalSlideReader):
    """Read KFBio-style IBL SQLite slides in physical coordinates.

    IBL stores every acquired 1664x1392 field as sixteen full-resolution JPEG
    tiles (layer 0) and one 4x downsampled JPEG (layer 1).  ``tbl_img_info``
    provides the field's level-0 ``nX/nY`` position.  A separate layer-2
    shrink table supplies a compact 16x overview used for the tissue mask.
    """

    def __init__(self, path: Path):
        self.path = Path(path)
        uri = f"file:{self.path.as_posix()}?mode=ro&immutable=1"
        self.connection = sqlite3.connect(uri, uri=True, check_same_thread=False)
        base = self.connection.execute(
            "SELECT total_img_width,total_img_height,img_width,img_height,"
            "tile_width,tile_height,ratio_step,pixel_size FROM tbl_base_info"
        ).fetchone()
        if base is None:
            raise ValueError(f"IBL base metadata is missing in {path}")
        (
            width,
            height,
            self._field_width,
            self._field_height,
            self._tile_width,
            self._tile_height,
            self._ratio_step,
            pixel_size,
        ) = base
        # The vendor stores millimetres per pixel (about 0.00025 at 40x).
        mpp = float(pixel_size) * 1000.0
        if not 0.05 <= mpp <= 2.0:
            raise ValueError(f"Implausible IBL physical pixel size {pixel_size!r} in {path}")
        records = self.connection.execute(
            "SELECT id,nX,nY FROM tbl_img_info ORDER BY id"
        ).fetchall()
        if not records:
            raise ValueError(f"IBL image-position table is empty in {path}")
        self._records = [
            (
                int(image_id),
                int(x0),
                int(y0),
                int(x0) + int(self._field_width),
                int(y0) + int(self._field_height),
            )
            for image_id, x0, y0 in records
        ]
        self._boxes = np.asarray([record[1:] for record in self._records], dtype=np.int64)
        self._spatial_bucket_px = max(int(self._field_width), int(self._field_height))
        buckets: Dict[Tuple[int, int], List[int]] = {}
        for index, (_, x0, y0, x1, y1) in enumerate(self._records):
            for bucket_y in range(y0 // self._spatial_bucket_px, (y1 - 1) // self._spatial_bucket_px + 1):
                for bucket_x in range(x0 // self._spatial_bucket_px, (x1 - 1) // self._spatial_bucket_px + 1):
                    buckets.setdefault((bucket_x, bucket_y), []).append(index)
        self._spatial_buckets = {
            key: np.asarray(sorted(set(indices)), dtype=np.int64)
            for key, indices in buckets.items()
        }
        self._field_cache: OrderedDict[Tuple[int, int], np.ndarray] = OrderedDict()
        self._field_cache_max = 320
        self.geometry = SlideGeometry(
            path=str(path),
            reader="ibl_sqlite",
            native_mpp_x=mpp,
            native_mpp_y=mpp,
            native_origin_x_px=0,
            native_origin_y_px=0,
            width_px=int(width),
            height_px=int(height),
        )

    @staticmethod
    def _decode_jpeg(blob: bytes) -> np.ndarray:
        bgr = cv2.imdecode(np.frombuffer(blob, dtype=np.uint8), cv2.IMREAD_COLOR)
        if bgr is None:
            raise ValueError("IBL JPEG tile could not be decoded")
        return cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)

    def _cached_field(self, image_id: int, layer: int) -> np.ndarray | None:
        key = (int(layer), int(image_id))
        if key in self._field_cache:
            value = self._field_cache.pop(key)
            self._field_cache[key] = value
            return value
        if layer == 1:
            row = self.connection.execute(
                "SELECT data FROM tbl_tile_info WHERE id=? AND layer=1 LIMIT 1",
                (int(image_id),),
            ).fetchone()
            if row is None:
                return None
            value = self._decode_jpeg(row[0])
        else:
            value = np.full(
                (int(self._field_height), int(self._field_width), 3),
                255,
                dtype=np.uint8,
            )
            found = False
            for col, row_index, blob in self.connection.execute(
                "SELECT col,row,data FROM tbl_tile_info "
                "WHERE id=? AND layer=0 ORDER BY row,col",
                (int(image_id),),
            ):
                tile = self._decode_jpeg(blob)
                x0 = int(col) * int(self._tile_width)
                y0 = int(row_index) * int(self._tile_height)
                x1 = min(value.shape[1], x0 + tile.shape[1])
                y1 = min(value.shape[0], y0 + tile.shape[0])
                if x1 > x0 and y1 > y0:
                    value[y0:y1, x0:x1] = tile[: y1 - y0, : x1 - x0]
                    found = True
            if not found:
                return None
        self._field_cache[key] = value
        while len(self._field_cache) > self._field_cache_max:
            self._field_cache.popitem(last=False)
        return value

    def _query_record_indices(self, x0: int, y0: int, x1: int, y1: int) -> np.ndarray:
        candidates: List[np.ndarray] = []
        for bucket_y in range(y0 // self._spatial_bucket_px, (y1 - 1) // self._spatial_bucket_px + 1):
            for bucket_x in range(x0 // self._spatial_bucket_px, (x1 - 1) // self._spatial_bucket_px + 1):
                found = self._spatial_buckets.get((bucket_x, bucket_y))
                if found is not None:
                    candidates.append(found)
        if not candidates:
            return np.empty(0, dtype=np.int64)
        indices = np.unique(np.concatenate(candidates))
        boxes = self._boxes[indices]
        keep = (
            (boxes[:, 2] > x0)
            & (boxes[:, 0] < x1)
            & (boxes[:, 3] > y0)
            & (boxes[:, 1] < y1)
        )
        return indices[keep]

    def _prefetch_fields(self, image_ids: Iterable[int], layer: int) -> None:
        missing = sorted(
            {
                int(image_id)
                for image_id in image_ids
                if (int(layer), int(image_id)) not in self._field_cache
            }
        )
        if not missing:
            return
        # SQLite stores the layer-1 JPEGs by image id.  Reading one ordered SQL
        # batch avoids hundreds of random BLOB lookups on large IBL files.
        for start in range(0, len(missing), 800):
            block = missing[start : start + 800]
            placeholders = ",".join("?" for _ in block)
            query = (
                "SELECT id,data FROM tbl_tile_info WHERE layer=? AND id IN ("
                + placeholders
                + ") ORDER BY id"
            )
            for image_id, blob in self.connection.execute(query, [int(layer), *block]):
                key = (int(layer), int(image_id))
                self._field_cache[key] = self._decode_jpeg(blob)
        while len(self._field_cache) > self._field_cache_max:
            self._field_cache.popitem(last=False)

    def read_center_patches(
        self,
        centers_um: np.ndarray,
        target_mpp: float,
        size: int,
        num_workers: int = 1,
        batch_size: int = 32,
    ) -> np.ndarray:
        del num_workers, batch_size
        layer = 1 if target_mpp >= 0.8 * self.geometry.native_mpp_x * self._ratio_step else 0
        fov = float(target_mpp) * int(size)
        required: set[int] = set()
        for cx, cy in np.asarray(centers_um, dtype=float):
            x0 = int(math.floor(max(0.0, cx - fov / 2.0) / self.geometry.native_mpp_x))
            y0 = int(math.floor(max(0.0, cy - fov / 2.0) / self.geometry.native_mpp_y))
            x1 = int(math.ceil(min(self.geometry.width_um, cx + fov / 2.0) / self.geometry.native_mpp_x))
            y1 = int(math.ceil(min(self.geometry.height_um, cy + fov / 2.0) / self.geometry.native_mpp_y))
            for index in self._query_record_indices(x0, y0, x1, y1):
                required.add(int(self._records[int(index)][0]))
        self._prefetch_fields(required, layer)
        result = []
        for cx, cy in np.asarray(centers_um, dtype=float):
            result.append(
                self.read_bbox_um(
                    (cx - fov / 2.0, cy - fov / 2.0, cx + fov / 2.0, cy + fov / 2.0),
                    target_mpp,
                )
            )
        return np.stack(result)

    def read_bbox_um(
        self, bbox_um: Tuple[float, float, float, float], target_mpp: float
    ) -> np.ndarray:
        x0_um, y0_um, x1_um, y1_um = map(float, bbox_um)
        x0_um, y0_um = max(0.0, x0_um), max(0.0, y0_um)
        x1_um = min(self.geometry.width_um, x1_um)
        y1_um = min(self.geometry.height_um, y1_um)
        out_w = max(1, int(round((x1_um - x0_um) / target_mpp)))
        out_h = max(1, int(round((y1_um - y0_um) / target_mpp)))
        x0 = int(math.floor(x0_um / self.geometry.native_mpp_x))
        y0 = int(math.floor(y0_um / self.geometry.native_mpp_y))
        x1 = int(math.ceil(x1_um / self.geometry.native_mpp_x))
        y1 = int(math.ceil(y1_um / self.geometry.native_mpp_y))
        # The patch pipeline requests 1.0 um/px. The vendor layer-1 image is 4*mpp
        # (typically 1.01 um/px) and is therefore the closest native level.
        layer = 1 if target_mpp >= 0.8 * self.geometry.native_mpp_x * self._ratio_step else 0
        downsample = int(self._ratio_step) ** int(layer)
        canvas = np.full((out_h, out_w, 3), 255, dtype=np.uint8)
        coverage = np.zeros((out_h, out_w), dtype=np.uint8)
        for index in self._query_record_indices(x0, y0, x1, y1):
            image_id, bx0, by0, bx1, by1 = self._records[int(index)]
            ix0, iy0 = max(x0, bx0), max(y0, by0)
            ix1, iy1 = min(x1, bx1), min(y1, by1)
            if ix1 <= ix0 or iy1 <= iy0:
                continue
            field = self._cached_field(image_id, layer)
            if field is None:
                continue
            sx0 = max(0, int(math.floor((ix0 - bx0) / downsample)))
            sy0 = max(0, int(math.floor((iy0 - by0) / downsample)))
            sx1 = min(field.shape[1], int(math.ceil((ix1 - bx0) / downsample)))
            sy1 = min(field.shape[0], int(math.ceil((iy1 - by0) / downsample)))
            if sx1 <= sx0 or sy1 <= sy0:
                continue
            dx0 = max(0, int(math.floor((ix0 * self.geometry.native_mpp_x - x0_um) / target_mpp)))
            dy0 = max(0, int(math.floor((iy0 * self.geometry.native_mpp_y - y0_um) / target_mpp)))
            dx1 = min(out_w, int(math.ceil((ix1 * self.geometry.native_mpp_x - x0_um) / target_mpp)))
            dy1 = min(out_h, int(math.ceil((iy1 * self.geometry.native_mpp_y - y0_um) / target_mpp)))
            if dx1 <= dx0 or dy1 <= dy0:
                continue
            region = cv2.resize(
                field[sy0:sy1, sx0:sx1],
                (dx1 - dx0, dy1 - dy0),
                interpolation=cv2.INTER_AREA if layer == 0 else cv2.INTER_CUBIC,
            )
            target_cov = coverage[dy0:dy1, dx0:dx1]
            write = target_cov == 0
            target = canvas[dy0:dy1, dx0:dx1]
            target[write] = region[write]
            target_cov[write] = 1
        return canvas

    def coarse_tissue_presence(self, target_mpp: float = 8.0) -> np.ndarray:
        rows = self.connection.execute(
            "SELECT layerNo,x,y,data FROM tbl_shrink_info ORDER BY y,x"
        ).fetchall()
        if not rows:
            return super().coarse_tissue_presence(target_mpp=target_mpp)
        layer_no = int(rows[0][0])
        downsample = int(self._ratio_step) ** layer_no
        width = max(1, int(math.ceil(self.geometry.width_px / downsample)))
        height = max(1, int(math.ceil(self.geometry.height_px / downsample)))
        overview = np.full((height, width, 3), 255, dtype=np.uint8)
        for _, x, y, blob in rows:
            tile = self._decode_jpeg(blob)
            tx0 = int(round(int(x) / downsample))
            ty0 = int(round(int(y) / downsample))
            tx1 = min(width, tx0 + tile.shape[1])
            ty1 = min(height, ty0 + tile.shape[0])
            if tx1 > tx0 and ty1 > ty0:
                overview[ty0:ty1, tx0:tx1] = tile[: ty1 - ty0, : tx1 - tx0]
        out_w = max(1, int(math.ceil(self.geometry.width_um / target_mpp)))
        out_h = max(1, int(math.ceil(self.geometry.height_um / target_mpp)))
        if overview.shape[:2] != (out_h, out_w):
            overview = cv2.resize(overview, (out_w, out_h), interpolation=cv2.INTER_AREA)
        array = overview.astype(np.float32)
        brightness = array.mean(axis=2)
        chroma = array.max(axis=2) - array.min(axis=2)
        return (brightness < 238.0) & (chroma > 6.0)

    def close(self) -> None:
        self._field_cache.clear()
        self.connection.close()


def _czi_scaling_um(metadata_xml: str) -> Tuple[float, float]:
    root = ET.fromstring(metadata_xml)
    values: Dict[str, float] = {}
    for distance in root.iter():
        if not distance.tag.endswith("Distance"):
            continue
        axis = distance.attrib.get("Id")
        value = None
        for child in distance.iter():
            if child.tag.endswith("Value") and child.text:
                value = float(child.text)
                break
        if axis in {"X", "Y"} and value is not None:
            values[axis] = value * 1e6
    if "X" not in values or "Y" not in values:
        raise ValueError("CZI metadata does not contain X/Y physical scaling")
    return values["X"], values["Y"]


class CziPhysicalReader(PhysicalSlideReader):
    def __init__(self, path: Path):
        import czifile
        from czifile import czifile as czi_module

        self.path = Path(path)
        self.czi = czifile.CziFile(str(path))
        mpp_x, mpp_y = _czi_scaling_um(self.czi.metadata())
        self._directory_recovery = "primary_directory"
        try:
            directory = list(self.czi.filtered_subblock_directory)
        except (UnicodeDecodeError, czi_module.SegmentNotFoundError):
            # A small number of SR386 files have an unreadable terminal directory
            # even though their image subblocks are intact.  Sequential segment
            # scanning follows the fallback already used internally by czifile,
            # but also handles a non-zero invalid directory pointer.
            scanned = [
                segment.directory_entry
                for segment in self.czi.segments(czi_module.SubBlockSegment.SID)
            ]
            mosaics = [entry for entry in scanned if entry.mosaic_index is not None]
            directory = sorted(mosaics, key=lambda entry: entry.mosaic_index) if mosaics else scanned
            self._directory_recovery = "sequential_subblock_scan"
        blocks = [sb for sb in directory if str(sb.pyramid_type) == "0"]
        if not blocks:
            raise ValueError(f"No full-resolution CZI subblocks in {path}")
        records = []
        for sb in blocks:
            dimensions = {entry.dimension: entry for entry in sb.dimension_entries}
            if "X" not in dimensions or "Y" not in dimensions:
                continue
            x_entry, y_entry = dimensions["X"], dimensions["Y"]
            x0, y0 = int(x_entry.start), int(y_entry.start)
            w, h = int(x_entry.size), int(y_entry.size)
            if w <= 0 or h <= 0:
                continue
            records.append((x0, y0, x0 + w, y0 + h, sb))
        if not records:
            raise ValueError(f"No readable X/Y full-resolution CZI subblocks in {path}")
        self._records = records
        self._boxes = np.asarray([record[:4] for record in records], dtype=np.int64)
        self._block_cache: OrderedDict[int, np.ndarray] = OrderedDict()
        self._unreadable_blocks: set[int] = set()
        # Row-major overlapping windows revisit the same CZI tiles on adjacent
        # grid rows.  Keep up to roughly 3 GiB of decoded tiles so a wide slide
        # does not repeatedly decompress every tile on the next row.
        median_block_bytes = int(
            np.median((self._boxes[:, 2] - self._boxes[:, 0]) * (self._boxes[:, 3] - self._boxes[:, 1]))
            * 3
        )
        self._block_cache_max = max(
            64,
            min(512, int((3 * 1024**3) // max(median_block_bytes, 1))),
        )
        # Avoid testing every one of tens of thousands of directory entries for
        # each 112-um patch.  This exact uniform-grid index only narrows the
        # candidate set; the final overlap predicate remains unchanged.
        self._spatial_bucket_px = max(
            1024,
            int(
                np.median(
                    np.maximum(
                        self._boxes[:, 2] - self._boxes[:, 0],
                        self._boxes[:, 3] - self._boxes[:, 1],
                    )
                )
            ),
        )
        spatial_buckets: Dict[Tuple[int, int], List[int]] = {}
        for index, (x0, y0, x1, y1) in enumerate(self._boxes):
            for bucket_y in range(y0 // self._spatial_bucket_px, (y1 - 1) // self._spatial_bucket_px + 1):
                for bucket_x in range(x0 // self._spatial_bucket_px, (x1 - 1) // self._spatial_bucket_px + 1):
                    spatial_buckets.setdefault((bucket_x, bucket_y), []).append(index)
        self._spatial_buckets = {
            key: np.asarray(sorted(set(indices)), dtype=np.int64)
            for key, indices in spatial_buckets.items()
        }
        min_x = int(self._boxes[:, 0].min())
        min_y = int(self._boxes[:, 1].min())
        max_x = int(self._boxes[:, 2].max())
        max_y = int(self._boxes[:, 3].max())
        self.geometry = SlideGeometry(
            path=str(path),
            reader=(
                "czifile_subblocks"
                if self._directory_recovery == "primary_directory"
                else "czifile_subblocks_recovered_directory"
            ),
            native_mpp_x=mpp_x,
            native_mpp_y=mpp_y,
            native_origin_x_px=min_x,
            native_origin_y_px=min_y,
            width_px=max_x - min_x,
            height_px=max_y - min_y,
        )

    @staticmethod
    def _rgb_from_subblock(sb) -> np.ndarray:
        array = np.asarray(sb.data_segment().data())
        while array.ndim > 3:
            array = array[0]
        if array.shape[-1] > 3:
            array = array[..., :3]
        return array.astype(np.uint8, copy=False)

    def _cached_rgb(self, index: int) -> np.ndarray | None:
        if index in self._unreadable_blocks:
            return None
        if index in self._block_cache:
            value = self._block_cache.pop(index)
            self._block_cache[index] = value
            return value
        try:
            value = self._rgb_from_subblock(self._records[index][4])
        except Exception as error:
            if error.__class__.__name__ != "SegmentNotFoundError" and "ZISRAW segment" not in str(error):
                raise
            self._unreadable_blocks.add(index)
            return None
        self._block_cache[index] = value
        while len(self._block_cache) > self._block_cache_max:
            self._block_cache.popitem(last=False)
        return value

    def _query_record_indices(self, x0: int, y0: int, x1: int, y1: int) -> np.ndarray:
        candidates: List[np.ndarray] = []
        for bucket_y in range(y0 // self._spatial_bucket_px, (y1 - 1) // self._spatial_bucket_px + 1):
            for bucket_x in range(x0 // self._spatial_bucket_px, (x1 - 1) // self._spatial_bucket_px + 1):
                found = self._spatial_buckets.get((bucket_x, bucket_y))
                if found is not None:
                    candidates.append(found)
        if not candidates:
            return np.empty(0, dtype=np.int64)
        candidate_index = np.unique(np.concatenate(candidates))
        boxes = self._boxes[candidate_index]
        keep = (
            (boxes[:, 2] > x0)
            & (boxes[:, 0] < x1)
            & (boxes[:, 3] > y0)
            & (boxes[:, 1] < y1)
        )
        return candidate_index[keep]

    def read_bbox_um(
        self, bbox_um: Tuple[float, float, float, float], target_mpp: float
    ) -> np.ndarray:
        x0_um, y0_um, x1_um, y1_um = bbox_um
        x0_um = max(0.0, x0_um)
        y0_um = max(0.0, y0_um)
        x1_um = min(self.geometry.width_um, x1_um)
        y1_um = min(self.geometry.height_um, y1_um)
        native_x0 = int(math.floor(x0_um / self.geometry.native_mpp_x))
        native_y0 = int(math.floor(y0_um / self.geometry.native_mpp_y))
        native_x1 = int(math.ceil(x1_um / self.geometry.native_mpp_x))
        native_y1 = int(math.ceil(y1_um / self.geometry.native_mpp_y))
        absolute_x0 = native_x0 + self.geometry.native_origin_x_px
        absolute_y0 = native_y0 + self.geometry.native_origin_y_px
        absolute_x1 = native_x1 + self.geometry.native_origin_x_px
        absolute_y1 = native_y1 + self.geometry.native_origin_y_px

        width = max(1, native_x1 - native_x0)
        height = max(1, native_y1 - native_y0)
        out_w = max(1, int(round((x1_um - x0_um) / target_mpp)))
        out_h = max(1, int(round((y1_um - y0_um) / target_mpp)))
        # A whole CZI can contain billions of native pixels but only a small
        # 8-um tissue-mask raster.  Assemble directly at target resolution in
        # that case instead of allocating the full native mosaic first.
        if width * height > 50_000_000 and target_mpp >= 4.0 * max(
            self.geometry.native_mpp_x, self.geometry.native_mpp_y
        ):
            canvas = np.full((out_h, out_w, 3), 255, dtype=np.uint8)
            coverage = np.zeros((out_h, out_w), dtype=np.uint8)
            hit = self._query_record_indices(absolute_x0, absolute_y0, absolute_x1, absolute_y1)
            origin_x = self.geometry.native_origin_x_px
            origin_y = self.geometry.native_origin_y_px
            for index in hit:
                bx0, by0, bx1, by1, _ = self._records[int(index)]
                ix0, iy0 = max(bx0, absolute_x0), max(by0, absolute_y0)
                ix1, iy1 = min(bx1, absolute_x1), min(by1, absolute_y1)
                if ix1 <= ix0 or iy1 <= iy0:
                    continue
                block = self._cached_rgb(int(index))
                if block is None:
                    continue
                region = block[iy0 - by0 : iy1 - by0, ix0 - bx0 : ix1 - bx0]
                physical_x0 = (ix0 - origin_x) * self.geometry.native_mpp_x
                physical_y0 = (iy0 - origin_y) * self.geometry.native_mpp_y
                physical_x1 = (ix1 - origin_x) * self.geometry.native_mpp_x
                physical_y1 = (iy1 - origin_y) * self.geometry.native_mpp_y
                dst_x0 = max(0, int(math.floor((physical_x0 - x0_um) / target_mpp)))
                dst_y0 = max(0, int(math.floor((physical_y0 - y0_um) / target_mpp)))
                dst_x1 = min(out_w, int(math.ceil((physical_x1 - x0_um) / target_mpp)))
                dst_y1 = min(out_h, int(math.ceil((physical_y1 - y0_um) / target_mpp)))
                if dst_x1 <= dst_x0 or dst_y1 <= dst_y0:
                    continue
                reduced = cv2.resize(
                    region,
                    (dst_x1 - dst_x0, dst_y1 - dst_y0),
                    interpolation=cv2.INTER_AREA,
                )
                target_cov = coverage[dst_y0:dst_y1, dst_x0:dst_x1]
                write = target_cov == 0
                target = canvas[dst_y0:dst_y1, dst_x0:dst_x1]
                target[write] = reduced[write]
                target_cov[write] = 1
            return canvas

        canvas = np.full((height, width, 3), 255, dtype=np.uint8)
        coverage = np.zeros((height, width), dtype=np.uint8)
        hit = self._query_record_indices(absolute_x0, absolute_y0, absolute_x1, absolute_y1)
        for index in hit:
            bx0, by0, bx1, by1, sb = self._records[int(index)]
            ix0, iy0 = max(bx0, absolute_x0), max(by0, absolute_y0)
            ix1, iy1 = min(bx1, absolute_x1), min(by1, absolute_y1)
            if ix1 <= ix0 or iy1 <= iy0:
                continue
            block = self._cached_rgb(int(index))
            if block is None:
                continue
            src_x0, src_y0 = ix0 - bx0, iy0 - by0
            src_x1, src_y1 = ix1 - bx0, iy1 - by0
            dst_x0, dst_y0 = ix0 - absolute_x0, iy0 - absolute_y0
            dst_x1, dst_y1 = ix1 - absolute_x0, iy1 - absolute_y0
            region = block[src_y0:src_y1, src_x0:src_x1]
            # Some CZI scenes overlap.  The first valid full-resolution block is
            # retained deterministically instead of averaging stain colours.
            target_cov = coverage[dst_y0:dst_y1, dst_x0:dst_x1]
            write = target_cov == 0
            target = canvas[dst_y0:dst_y1, dst_x0:dst_x1]
            target[write] = region[write]
            target_cov[write] = 1

        if canvas.shape[:2] != (out_h, out_w):
            interpolation = cv2.INTER_AREA if out_w <= width and out_h <= height else cv2.INTER_CUBIC
            canvas = cv2.resize(canvas, (out_w, out_h), interpolation=interpolation)
        return canvas

    def coarse_tissue_presence(self, target_mpp: float = 8.0) -> np.ndarray:
        """Use the sparse CZI tile directory to avoid decoding a whole slide.

        This is only a candidate-presence mask.  Every UNI tile is subsequently
        checked with the RGB tissue fraction, so blank acquired subblocks do not
        enter the model merely because their directory entry exists.
        """

        width = max(1, int(math.ceil(self.geometry.width_um / target_mpp)))
        height = max(1, int(math.ceil(self.geometry.height_um / target_mpp)))
        mask = np.zeros((height, width), dtype=np.uint8)
        origin_x = self.geometry.native_origin_x_px
        origin_y = self.geometry.native_origin_y_px
        for x0, y0, x1, y1, _ in self._records:
            tx0 = max(0, int(math.floor((x0 - origin_x) * self.geometry.native_mpp_x / target_mpp)))
            ty0 = max(0, int(math.floor((y0 - origin_y) * self.geometry.native_mpp_y / target_mpp)))
            tx1 = min(width, int(math.ceil((x1 - origin_x) * self.geometry.native_mpp_x / target_mpp)))
            ty1 = min(height, int(math.ceil((y1 - origin_y) * self.geometry.native_mpp_y / target_mpp)))
            if tx1 > tx0 and ty1 > ty0:
                mask[ty0:ty1, tx0:tx1] = 1
        return mask.astype(bool)

    def close(self) -> None:
        self._block_cache.clear()
        self.czi.close()


def open_physical_slide(
    path: Path, mpp_x: float | None = None, mpp_y: float | None = None
) -> PhysicalSlideReader:
    suffix = path.suffix.lower()
    if suffix == ".ibl":
        return IblPhysicalReader(path)
    if suffix == ".czi":
        return CziPhysicalReader(path)
    return OpenSlidePhysicalReader(path, mpp_x=mpp_x, mpp_y=mpp_y)


def tissue_fraction(rgb: np.ndarray) -> float:
    array = rgb.astype(np.float32)
    brightness = array.mean(axis=2)
    chroma = array.max(axis=2) - array.min(axis=2)
    tissue = (brightness < 232.0) & (chroma > 8.0)
    return float(tissue.mean())
