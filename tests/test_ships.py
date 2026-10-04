import pytest

from underway import transfer
from underway.ships import SHIPS
from underway.source import Source


def test_ships_table_has_four_ships():
    assert sorted(SHIPS) == ["armstrong", "discovery", "revelle", "sikuliaq"]


@pytest.mark.parametrize("ship", sorted(SHIPS))
def test_ship_source_names_unique(ship):
    names = [s.name for s in SHIPS[ship].sources]
    assert len(names) == len(set(names))


@pytest.mark.parametrize("ship", sorted(SHIPS))
def test_ship_source_drives_have_a_server(ship):
    drives = {s.drive for s in SHIPS[ship].sources}
    assert drives <= set(SHIPS[ship].servers)


@pytest.mark.parametrize("ship", sorted(SHIPS))
def test_ship_source_transfer_methods_exist(ship):
    assert {s.transfer for s in SHIPS[ship].sources} <= set(transfer.TRANSFER)


@pytest.mark.parametrize("ship", sorted(SHIPS))
def test_ship_met_sources_have_parsers(ship):
    assert all(s.parser is not None for s in SHIPS[ship].sources if s.met)


def test_source_defaults():
    src = Source("x", drive="d", remote="r")
    assert (
        src.pattern,
        src.transfer,
        src.parser,
        src.cache,
        src.met,
        src.readonly,
        src.local,
    ) == (
        "*",
        "rsync",
        None,
        True,
        False,
        False,
        None,
    )
