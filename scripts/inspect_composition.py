#!/usr/bin/env python3
import json
import sys
from pathlib import Path

report = json.loads(Path(sys.argv[1]).read_text())
print(json.dumps(report.get("complex_composition", {}), indent=2))
