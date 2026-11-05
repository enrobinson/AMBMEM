def test_example_disables_reflection():
    import yaml
    config = yaml.safe_load(open("examples/membrane_only_custom_sterol.yaml"))
    species = config["custom_sterols"]["species"][0]
    assert species["orientation"]["allow_reflection"] is False
