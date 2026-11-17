"""Persistence codec for supercell results.

The envelope mirrors the primitive result codec (version plus result), but the
result carries the integer matrix, the derived canonical cell, the origin map
and the full ledger.  Decoding only validates and rebuilds value objects; it
never enumerates replicas or reconstructs the supercell from the source.
"""
from dataclasses import dataclass
from fractions import Fraction
from .exact_numbers import (
    decode_rational, encode_rational, exact_keys, json_dump, json_load,
    read_identity, read_integer)
from .bonds import read_integer_offset
from .bond_geometry import BondGeometry
from .input_model import cell_to_mapping
from .intake import accept_cell
from .result_codec import (
    SUPERCELL_CALCULATION, read_source, source_to_mapping)
from .supercell_matrix import read_supercell_matrix

SUPERCELL_RESULT_VERSION = 1
UNITS = {"length": "nm", "length_squared": "nm2", "volume": "nm3"}


@dataclass(frozen=True)
class SupercellResult:
    source: tuple
    matrix: tuple
    cell: object
    origin_sites: tuple
    origin_bonds: tuple
    volume: Fraction
    positions: tuple
    bonds: tuple


def supercell_to_mapping(result):
    return {
        "source": source_to_mapping(result.source),
        "units": dict(UNITS),
        "matrix": [list(row) for row in result.matrix],
        "cell": cell_to_mapping(result.cell),
        "origins": {
            "sites": [{"id": identity, "source_id": source_id,
                       "translation": list(translation)}
                      for identity, source_id, translation in result.origin_sites],
            "bonds": [{"id": identity, "source_id": source_id}
                      for identity, source_id in result.origin_bonds],
        },
        "volume": result.volume,
        "positions": [{"id": identity, "cartesian": list(position)}
                      for identity, position in result.positions],
        "bonds": [{"id": row.id, "u": row.u, "v": row.v,
                   "offset": list(row.offset),
                   "displacement": list(row.displacement),
                   "length_squared": row.length_squared}
                  for row in result.bonds],
    }


def encode_supercell_result(result):
    value = supercell_to_mapping(result)
    value["volume"] = encode_rational(value["volume"])
    for site in value["positions"]:
        site["cartesian"] = [encode_rational(x) for x in site["cartesian"]]
    for bond in value["bonds"]:
        bond["displacement"] = [encode_rational(x) for x in bond["displacement"]]
        bond["length_squared"] = encode_rational(bond["length_squared"])
    value["matrix"] = [[read_integer(x, "matrix entry") for x in row]
                       for row in value["matrix"]]
    value["origins"]["sites"] = [
        {"id": item["id"], "source_id": item["source_id"],
         "translation": [read_integer(x, "translation") for x in item["translation"]]}
        for item in value["origins"]["sites"]]
    text = json_dump({"version": SUPERCELL_RESULT_VERSION, "result": value})
    decode_supercell_result(text)
    return text


def _vector(value):
    if type(value) is not list or len(value) != 3:
        raise ValueError("incomplete supercell vector")
    return tuple(decode_rational(x) for x in value)


def _origins_sites(value, site_ids):
    if type(value) is not list:
        raise ValueError("incomplete site origins")
    rows, seen = [], set()
    for item in value:
        exact_keys(item, ("id", "source_id", "translation"))
        identity = read_identity(item["id"])
        source_id = read_identity(item["source_id"])
        translation = read_integer_offset(item["translation"])
        if identity in seen or identity not in site_ids:
            raise ValueError("invalid site origin identity")
        seen.add(identity)
        rows.append((identity, source_id, translation))
    if seen != site_ids:
        raise ValueError("site origins do not cover the derived cell")
    return tuple(sorted(rows))


def _origins_bonds(value, bond_ids):
    if type(value) is not list:
        raise ValueError("incomplete bond origins")
    rows, seen = [], set()
    for item in value:
        exact_keys(item, ("id", "source_id"))
        identity = read_identity(item["id"])
        source_id = read_identity(item["source_id"])
        if identity in seen or identity not in bond_ids:
            raise ValueError("invalid bond origin identity")
        seen.add(identity)
        rows.append((identity, source_id))
    if seen != bond_ids:
        raise ValueError("bond origins do not cover the derived ledger")
    return tuple(sorted(rows))


def decode_supercell_result(text):
    payload = json_load(text)
    exact_keys(payload, ("version", "result"))
    if type(payload["version"]) is not int or payload["version"] != 1:
        raise ValueError("unsupported supercell result format")
    value = exact_keys(payload["result"], (
        "source", "units", "matrix", "cell", "origins",
        "volume", "positions", "bonds"))
    source = read_source(value["source"], SUPERCELL_CALCULATION)
    if value["units"] != UNITS:
        raise ValueError("incorrect supercell result units")
    matrix = read_supercell_matrix(value["matrix"])
    cell, _ = accept_cell(value["cell"])
    if cell_to_mapping(cell) != value["cell"]:
        raise ValueError("derived supercell cell is not canonical")
    site_ids = {site.id for site in cell.sites}
    bond_ids = {bond.id for bond in cell.bonds}
    origins = exact_keys(value["origins"], ("sites", "bonds"))
    origin_sites = _origins_sites(origins["sites"], site_ids)
    origin_bonds = _origins_bonds(origins["bonds"], bond_ids)
    volume = decode_rational(value["volume"])
    if volume <= 0:
        raise ValueError("nonpositive supercell volume")
    if type(value["positions"]) is not list or type(value["bonds"]) is not list:
        raise ValueError("incomplete supercell ledger")
    positions, bonds = [], []
    position_ids, ledger_ids = set(), set()
    for site in value["positions"]:
        exact_keys(site, ("id", "cartesian"))
        identity = read_identity(site["id"])
        if identity in position_ids:
            raise ValueError("duplicate supercell site")
        position_ids.add(identity)
        positions.append((identity, _vector(site["cartesian"])))
    if position_ids != site_ids:
        raise ValueError("position coverage differs from derived cell")
    for row in value["bonds"]:
        exact_keys(row, ("id", "u", "v", "offset",
                         "displacement", "length_squared"))
        identity = read_identity(row["id"])
        u, v = read_identity(row["u"]), read_identity(row["v"])
        if identity in ledger_ids or u not in site_ids or v not in site_ids:
            raise ValueError("invalid supercell bond identity or endpoint")
        ledger_ids.add(identity)
        displacement = _vector(row["displacement"])
        length = decode_rational(row["length_squared"])
        if length != sum((x * x for x in displacement), Fraction(0)):
            raise ValueError("inconsistent supercell length")
        bonds.append(BondGeometry(identity, u, v,
                                  read_integer_offset(row["offset"]),
                                  displacement, length))
    if ledger_ids != bond_ids:
        raise ValueError("ledger coverage differs from derived cell")
    return SupercellResult(
        source, matrix, cell, origin_sites, origin_bonds, volume,
        tuple(sorted(positions)), tuple(sorted(bonds, key=lambda row: row.id)))
