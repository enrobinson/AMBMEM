from ambmem.cli import main


def test_import_cli():
    assert callable(main)
