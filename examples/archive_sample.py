"""Reproducible mini-example of the SpectrumArchive delivery workflow.

Runs entirely in a temporary directory.  Usage:

    PYTHONPATH=src python3 examples/archive_sample.py
"""

import json
import tempfile

from wavelength_calibration import SpectrumArchive

A_DATA = {
    "max_residual_nm": "0",
    "reference_lines": [
        {"id": "l0", "pixel": "0", "wavelength_nm": "400"},
        {"id": "l1", "pixel": "10", "wavelength_nm": "421"},
        {"id": "l2", "pixel": "20", "wavelength_nm": "444"},
        {"id": "l3", "pixel": "30", "wavelength_nm": "469"},
    ],
}


def show(title, value):
    print(f"== {title} ==")
    print(json.dumps(value, indent=2, ensure_ascii=False))
    print()


with tempfile.TemporaryDirectory() as root:
    arc = SpectrumArchive(root)

    calibration = arc.set_calibration("A", 0, A_DATA)
    show("calibration revision 1", calibration)

    arc.add_exposure(
        "EA",
        "A",
        pixel_edges=["0", "10", "20"],
        counts=["21", "23"],
        wavelength_edges=["400", "410", "421", "430", "444"],
    )
    product = arc.build("EA")
    show("product: exact counts [10, 11, 9, 14], dropped 0", product)

    release = arc.release("request-1", ["EA"])
    show("release number 1", release)

    # Correct the lineage: only a reference line id changes.  The exposure is
    # stale immediately even though every number would come out identical.
    renamed = json.loads(json.dumps(A_DATA))
    renamed["reference_lines"][2]["id"] = "l2-renamed"
    arc.set_calibration("A", 1, renamed)
    show("after correction: EA is stale", arc.inspect_exposure("EA")["state"])

    # The original request still replays the frozen delivery verbatim.
    replay = arc.release("request-1", ["EA"])
    show("replay returns the original release", replay["release_id"])

    # Rebuilding produces a new product identity at the new evidence version.
    rebuilt = arc.build("EA")
    show("rebuilt product at revision 2", {k: rebuilt[k] for k in ("product_id", "calibration_revision", "target_counts")})

    arc.close()
