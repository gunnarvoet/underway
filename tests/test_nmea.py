import logging

import numpy as np
import pandas as pd
import pytest

from underway.parsers import nmea


def _fix(talker="GP", zda_time="120000.00", gga_time="120000.00", ns="N", ew="W"):
    return [
        f"${talker}ZDA,{zda_time},26,11,2025,,",
        f"${talker}GGA,{gga_time},2118.600000,{ns},15752.200000,{ew},2,10,0.9,14.26,M",
        f"${talker}VTG,4.32,T,354.83,M,10.0,N,18.5,K,D",
        f"${talker}HDT,26.63,T",
    ]


def test_parse_southern_and_eastern_hemisphere_signs():
    ds = nmea.parse(pd.Series(_fix(ns="S", ew="E")))
    assert float(ds.lat[0]) == pytest.approx(-(21 + 18.6 / 60))
    assert float(ds.lon[0]) == pytest.approx(157 + 52.2 / 60)


def test_parse_other_talker_id_gives_position_and_heading():
    ds = nmea.parse(pd.Series(_fix(talker="IN")))
    assert float(ds.lat[0]) == pytest.approx(21 + 18.6 / 60)
    assert float(ds.heading[0]) == pytest.approx(26.63)


def test_parse_gga_time_without_decimals_matches_zda():
    ds = nmea.parse(pd.Series(_fix(gga_time="120000")))
    assert float(ds.lat[0]) == pytest.approx(21 + 18.6 / 60)


def test_parse_hour_out_of_range_fix_dropped():
    sentences = _fix(zda_time="996000.00", gga_time="996000.00") + _fix()
    ds = nmea.parse(pd.Series(sentences))
    assert ds.time.values.tolist() == [
        np.datetime64("2025-11-26T12:00:00", "ns").astype(int)
    ]


def test_parse_position_ahead_of_zda_warns_about_sentence_order(caplog):
    one, two = _fix(), _fix(zda_time="120001.00", gga_time="120001.00")
    sentences = [one[1], one[0], one[2], one[3], two[1], two[0], two[2], two[3]]
    with caplog.at_level(logging.WARNING, logger="underway"):
        nmea.parse(pd.Series(sentences), name="x")
    assert "sentence order" in caplog.text
