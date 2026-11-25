"""Count-conserving spectral product correction and traceable delivery.

This module adds a persistent archive on top of the exact quadratic fit in
:mod:`wavelength_calibration.service`.  It never re-implements the fit: every
calibration candidate is evaluated by calling :func:`service.process` with the
two endpoint reference pixels as the query pixels, and the resulting evidence
(coefficients, per-line residuals, maximum residual, acceptance status) is
stored verbatim with the calibration revision.

All arithmetic is exact rational arithmetic (:class:`fractions.Fraction`).
External numeric inputs are plain decimal strings of at most 32 characters;
outputs are reduced ``n`` or ``n/d`` strings inside plain JSON-serializable
objects.
"""

from __future__ import annotations

import copy
import json
import os
import re
from fractions import Fraction
from pathlib import Path

from .service import process

__all__ = ["SpectrumArchive"]

_DECIMAL = re.compile(r"[+-]?(\d+(\.\d*)?|\.\d+)\Z")
_MAX_N_BINS = 2000


def _strict_num(value: object, name: str) -> Fraction:
    """Parse a plain decimal string; reject numbers, NaN/inf, exponents, fractions."""
    if not isinstance(value, str) or not value or len(value) > 32 or not _DECIMAL.match(value):
        raise ValueError(f"{name} must be a plain decimal string of at most 32 characters")
    return Fraction(value)


def _text(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def _require_id(value: object, name: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a nonempty string")
    return value


def _require_revision(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


class SpectrumArchive:
    """Persistent archive of calibration lineages, raw exposures, products and releases.

    ``root`` is a dedicated directory owned by the archive.  Opening an empty
    directory creates a new archive; reopening the same directory continues
    where the previous close left off.  State is a single JSON document
    replaced atomically on every committed change.
    """

    FORMAT = 1
    _STATE_FILE = "state.json"

    def __init__(self, root, *, fault_hook=None):
        if fault_hook is not None and not callable(fault_hook):
            raise ValueError("fault_hook must be callable")
        self._hook = fault_hook
        self._closed = False
        self._root = Path(root)
        self._root.mkdir(parents=True, exist_ok=True)
        self._path = self._root / self._STATE_FILE
        if self._path.exists():
            with open(self._path, "r", encoding="utf-8") as fh:
                state = json.load(fh)
            if not isinstance(state, dict) or state.get("format") != self.FORMAT:
                raise ValueError("unsupported archive format")
            self._state = state
        else:
            # First open of an empty directory: initialize.  This is not a
            # committed business operation, so the fault hook is not invoked.
            self._state = {
                "format": self.FORMAT,
                "slots": {},
                "exposures": {},
                "releases": {},
                "requests": {},
                "next_release_number": 1,
                "next_product_seq": 1,
            }
            self._persist(self._state)

    # ------------------------------------------------------------------
    # persistence and commit protocol
    # ------------------------------------------------------------------

    def _persist(self, state: dict) -> None:
        tmp = self._root / (self._STATE_FILE + ".tmp")
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(state, fh, ensure_ascii=False, sort_keys=True)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, self._path)
        try:
            dfd = os.open(self._root, os.O_RDONLY | os.O_DIRECTORY)
        except (OSError, AttributeError):
            return
        try:
            os.fsync(dfd)
        finally:
            os.close(dfd)

    def _commit(self, new_state: dict, operation: str) -> None:
        """Make ``new_state`` durable, invoking the fault hook at the two boundaries."""
        if self._hook is not None:
            self._hook("before_commit", operation)
        self._persist(new_state)
        self._state = new_state
        if self._hook is not None:
            self._hook("after_commit", operation)

    def _check_open(self) -> None:
        if self._closed:
            raise ValueError("archive is closed")

    def close(self) -> None:
        self._closed = True

    def __enter__(self) -> "SpectrumArchive":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    # ------------------------------------------------------------------
    # calibration lineage
    # ------------------------------------------------------------------

    def set_calibration(self, slot_id: str, expected_revision: int, data: dict) -> dict:
        """Create or correct a calibration lineage.

        ``expected_revision`` must be 0 to create the lineage, or equal to the
        current integer revision to correct it.  On success the revision is
        incremented (starting from 1) and the new evidence record is returned.
        """
        self._check_open()
        _require_id(slot_id, "slot_id")
        _require_revision(expected_revision, "expected_revision")
        if not isinstance(data, dict):
            raise ValueError("data must be a mapping")
        lines = data.get("reference_lines")
        tolerance_raw = data.get("max_residual_nm")
        _strict_num(tolerance_raw, "max_residual_nm")
        if not isinstance(lines, list):
            raise ValueError("reference_lines must be a list")
        for line in lines:
            if not isinstance(line, dict):
                raise ValueError("each reference line must be a mapping")
            _require_id(line.get("id"), "reference line id")
            _strict_num(line.get("pixel"), "pixel")
            _strict_num(line.get("wavelength_nm"), "wavelength_nm")
        # Delegate the fit to the existing exact engine, querying the two
        # endpoint reference pixels so the full fit evidence is captured.
        query = [lines[0]["pixel"], lines[-1]["pixel"]] if lines else ["0"]
        result = process({
            "reference_lines": lines,
            "max_residual_nm": tolerance_raw,
            "query_pixels": query,
        })
        if result["status"] != "accepted":
            raise ValueError("calibration residuals exceed max_residual_nm")
        c0 = Fraction(result["coefficients"]["c0"])
        c1 = Fraction(result["coefficients"]["c1"])
        c2 = Fraction(result["coefficients"]["c2"])
        p_lo = Fraction(lines[0]["pixel"])
        p_hi = Fraction(lines[-1]["pixel"])
        # The derivative c1 + 2*c2*p is linear, so it is strictly positive on
        # the whole closed reference interval iff it is positive at both ends.
        if c1 + 2 * c2 * p_lo <= 0 or c1 + 2 * c2 * p_hi <= 0:
            raise ValueError("fit is not strictly increasing over the reference pixel range")
        # With a strictly positive derivative the mapped wavelength is minimal
        # at the left endpoint; it must be positive there.
        if c0 + c1 * p_lo + c2 * p_lo * p_lo <= 0:
            raise ValueError("fit maps the reference range to nonpositive wavelengths")
        slot = self._state["slots"].get(slot_id)
        current = slot["current_revision"] if slot is not None else 0
        if expected_revision != current:
            raise ValueError(
                f"expected_revision {expected_revision} does not match current revision {current}"
            )
        revision = current + 1
        record = {
            "slot_id": slot_id,
            "revision": revision,
            "input": {
                "reference_lines": copy.deepcopy(lines),
                "max_residual_nm": tolerance_raw,
            },
            "coefficients": result["coefficients"],
            "reference_residuals": result["reference_residuals"],
            "max_abs_residual_nm": result["max_abs_residual_nm"],
            "status": result["status"],
            "calibrated": result["calibrated"],
        }
        new_state = copy.deepcopy(self._state)
        entry = new_state["slots"].setdefault(slot_id, {"current_revision": 0, "revisions": {}})
        entry["revisions"][str(revision)] = record
        entry["current_revision"] = revision
        self._commit(new_state, "set_calibration")
        return copy.deepcopy(record)

    def get_calibration(self, slot_id: str, revision: int | None = None) -> dict:
        """Return a stored calibration record; never refits or rewrites history."""
        self._check_open()
        _require_id(slot_id, "slot_id")
        slot = self._state["slots"].get(slot_id)
        if slot is None:
            raise ValueError("unknown slot_id")
        if revision is None:
            revision = slot["current_revision"]
        else:
            _require_revision(revision, "revision")
        record = slot["revisions"].get(str(revision))
        if record is None:
            raise ValueError("unknown calibration revision")
        return copy.deepcopy(record)

    # ------------------------------------------------------------------
    # raw exposures
    # ------------------------------------------------------------------

    @staticmethod
    def _num_list(value: object, name: str) -> list[Fraction]:
        if not isinstance(value, list) or not value:
            raise ValueError(f"{name} must be a nonempty list of decimal strings")
        return [_strict_num(item, f"{name} item") for item in value]

    def add_exposure(
        self,
        exposure_id: str,
        slot_id: str,
        pixel_edges: list[str],
        counts: list[str],
        wavelength_edges: list[str],
    ) -> dict:
        """Register an immutable raw exposure bound to a calibration lineage."""
        self._check_open()
        _require_id(exposure_id, "exposure_id")
        _require_id(slot_id, "slot_id")
        if slot_id not in self._state["slots"]:
            raise ValueError("unknown slot_id")
        if exposure_id in self._state["exposures"]:
            raise ValueError("exposure_id is already registered")
        pixels = self._num_list(pixel_edges, "pixel_edges")
        counts_v = self._num_list(counts, "counts")
        waves = self._num_list(wavelength_edges, "wavelength_edges")
        if not 2 <= len(pixels) <= _MAX_N_BINS + 1:
            raise ValueError("pixel_edges must contain 2..2001 edges")
        if len(counts_v) != len(pixels) - 1:
            raise ValueError("counts must have exactly one fewer item than pixel_edges")
        if not 2 <= len(waves) <= _MAX_N_BINS + 1:
            raise ValueError("wavelength_edges must contain 2..2001 edges")
        for left, right in zip(pixels, pixels[1:]):
            if right <= left:
                raise ValueError("pixel_edges must be strictly increasing")
        for left, right in zip(waves, waves[1:]):
            if right <= left:
                raise ValueError("wavelength_edges must be strictly increasing")
        if pixels[0] < 0:
            raise ValueError("pixel edges must be nonnegative")
        if waves[0] <= 0:
            raise ValueError("wavelength edges must be positive")
        for count in counts_v:
            if count < 0:
                raise ValueError("counts must be nonnegative")
        record = {
            "exposure_id": exposure_id,
            "slot_id": slot_id,
            "pixel_edges": list(pixel_edges),
            "counts": list(counts),
            "wavelength_edges": list(wavelength_edges),
        }
        new_state = copy.deepcopy(self._state)
        new_state["exposures"][exposure_id] = {**copy.deepcopy(record), "products": {}}
        self._commit(new_state, "add_exposure")
        return record

    # ------------------------------------------------------------------
    # products
    # ------------------------------------------------------------------

    def _exposure_status(self, exposure: dict) -> tuple[str, dict | None]:
        slot = self._state["slots"][exposure["slot_id"]]
        products = exposure["products"]
        if not products:
            return "unbuilt", None
        last_rev = max(int(rev) for rev in products)
        last = products[str(last_rev)]
        state = "ready" if last_rev == slot["current_revision"] else "stale"
        return state, last

    def build(self, exposure_id: str) -> dict:
        """Build the wavelength-grid product for an exposure at the current calibration.

        Rebuilding an exposure that already has a product at the current
        calibration revision returns that product unchanged.  A successful
        calibration correction forces a new product identity on the next build.
        """
        self._check_open()
        _require_id(exposure_id, "exposure_id")
        exposure = self._state["exposures"].get(exposure_id)
        if exposure is None:
            raise ValueError("unknown exposure_id")
        slot = self._state["slots"][exposure["slot_id"]]
        revision = slot["current_revision"]
        existing = exposure["products"].get(str(revision))
        if existing is not None:
            return copy.deepcopy(existing)
        calibration = slot["revisions"][str(revision)]
        c0 = Fraction(calibration["coefficients"]["c0"])
        c1 = Fraction(calibration["coefficients"]["c1"])
        c2 = Fraction(calibration["coefficients"]["c2"])
        lines = calibration["input"]["reference_lines"]
        p_lo = Fraction(lines[0]["pixel"])
        p_hi = Fraction(lines[-1]["pixel"])
        pixels = [Fraction(s) for s in exposure["pixel_edges"]]
        if pixels[0] < p_lo or pixels[-1] > p_hi:
            raise ValueError("pixel edges are outside the calibrated reference pixel range")
        mapped = [c0 + c1 * p + c2 * p * p for p in pixels]
        targets = [Fraction(s) for s in exposure["wavelength_edges"]]
        if targets[0] < mapped[0] or targets[-1] > mapped[-1]:
            raise ValueError("wavelength edges are outside the mapped source wavelength range")
        counts = [Fraction(s) for s in exposure["counts"]]
        # Each source bin's integrated count is uniformly distributed over its
        # mapped wavelength interval; each target bin collects count * overlap
        # / source width, summed exactly.
        n_target = len(targets) - 1
        target_counts = [Fraction(0)] * n_target
        j = 0
        for i, count in enumerate(counts):
            lo, hi = mapped[i], mapped[i + 1]
            width = hi - lo
            while j < n_target and targets[j + 1] <= lo:
                j += 1
            jj = j
            while jj < n_target and targets[jj] < hi:
                overlap = min(hi, targets[jj + 1]) - max(lo, targets[jj])
                if overlap > 0:
                    target_counts[jj] += count * overlap / width
                jj += 1
        dropped = sum(counts, Fraction(0)) - sum(target_counts, Fraction(0))
        assert dropped >= 0, "conservation invariant violated"
        seq = self._state["next_product_seq"]
        product = {
            "product_id": f"P{seq:06d}",
            "exposure_id": exposure_id,
            "slot_id": exposure["slot_id"],
            "calibration_revision": revision,
            "pixel_edges": list(exposure["pixel_edges"]),
            "counts": list(exposure["counts"]),
            "wavelength_edges": list(exposure["wavelength_edges"]),
            "source_wavelength_edges": [_text(x) for x in mapped],
            "target_counts": [_text(x) for x in target_counts],
            "dropped_counts": _text(dropped),
            "calibration": copy.deepcopy(calibration),
        }
        new_state = copy.deepcopy(self._state)
        new_state["exposures"][exposure_id]["products"][str(revision)] = product
        new_state["next_product_seq"] = seq + 1
        self._commit(new_state, "build")
        return copy.deepcopy(product)

    def inspect_exposure(self, exposure_id: str) -> dict:
        """Return the raw exposure plus its delivery state and last product."""
        self._check_open()
        _require_id(exposure_id, "exposure_id")
        exposure = self._state["exposures"].get(exposure_id)
        if exposure is None:
            raise ValueError("unknown exposure_id")
        state, last = self._exposure_status(exposure)
        return {
            "exposure_id": exposure["exposure_id"],
            "slot_id": exposure["slot_id"],
            "pixel_edges": list(exposure["pixel_edges"]),
            "counts": list(exposure["counts"]),
            "wavelength_edges": list(exposure["wavelength_edges"]),
            "state": state,
            "last_product": copy.deepcopy(last),
        }

    # ------------------------------------------------------------------
    # releases
    # ------------------------------------------------------------------

    def release(self, request_id: str, exposure_ids: list[str]) -> dict:
        """Deliver the current ready products of a batch of exposures.

        ``request_id`` is an idempotent identity: replaying a succeeded request
        with the same exposure set returns the original delivery verbatim, even
        if the exposures have since become stale or been rebuilt.  A failed
        request consumes neither the request_id nor a release number.
        """
        self._check_open()
        _require_id(request_id, "request_id")
        if not isinstance(exposure_ids, list) or not exposure_ids:
            raise ValueError("exposure_ids must be a nonempty list")
        for exposure_id in exposure_ids:
            _require_id(exposure_id, "exposure_id")
        if len(set(exposure_ids)) != len(exposure_ids):
            raise ValueError("exposure_ids must not contain duplicates")
        exposures = self._state["exposures"]
        for exposure_id in exposure_ids:
            if exposure_id not in exposures:
                raise ValueError(f"unknown exposure_id {exposure_id!r}")
        key = sorted(exposure_ids)
        binding = self._state["requests"].get(request_id)
        if binding is not None:
            if binding["exposure_ids"] != key:
                raise ValueError("request_id is bound to a different exposure set")
            return copy.deepcopy(self._state["releases"][binding["release_id"]])
        products = []
        for exposure_id in key:
            state, last = self._exposure_status(exposures[exposure_id])
            if state != "ready":
                raise ValueError(f"exposure {exposure_id!r} is {state}, not ready")
            products.append(copy.deepcopy(last))
        number = self._state["next_release_number"]
        release = {
            "release_id": f"R{number:06d}",
            "release_number": number,
            "request_id": request_id,
            "products": products,
        }
        new_state = copy.deepcopy(self._state)
        new_state["releases"][release["release_id"]] = release
        new_state["requests"][request_id] = {
            "exposure_ids": key,
            "release_id": release["release_id"],
        }
        new_state["next_release_number"] = number + 1
        self._commit(new_state, "release")
        return copy.deepcopy(release)

    def get_release(self, release_id: str) -> dict:
        """Return a past delivery exactly as it was made."""
        self._check_open()
        _require_id(release_id, "release_id")
        release = self._state["releases"].get(release_id)
        if release is None:
            raise ValueError("unknown release_id")
        return copy.deepcopy(release)
