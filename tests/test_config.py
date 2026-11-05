from pathlib import Path
import yaml


def test_control_example_parses():
    config = yaml.safe_load(Path("examples/membrane_only_control.yaml").read_text())
    assert config["membrane_only"]["enabled"] is True
    assert config["membrane_only"]["target_per_leaflet"] == 200
