import json
from pathlib import Path
from wavelength_calibration import process

p=Path(__file__).parents[1]/"fixtures"/"input.json"
print(json.dumps(process(json.loads(p.read_text())),indent=2))

