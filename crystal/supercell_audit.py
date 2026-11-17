"""Independent exact audit of a persisted supercell result.

The audit consumes only the fixed source cell, the stored integer matrix and
the stored derived cell/origins/ledger.  It re-derives every expected number
with exact rational arithmetic (no floats, no tolerances) and rejects any
missing, extra or mismatched evidence.  Read and export call it, but it never
enumerates replicas, so a tampered body cannot be repaired by the audit.
"""
from fractions import Fraction
from .input_digest import verify_cell_digest
from .rational_matrix import (
    add3, determinant3, matmul3, matvec3, subtract3)
from .sites import index_sites
from .bonds import read_integer_offset
from .supercell_matrix import coset_key


def _check_site_origins(source_cell, matrix, result):
    source_sites = index_sites(source_cell.sites)
    derived_sites = index_sites(result.cell.sites)
    if len(result.origin_sites) != len(derived_sites):
        raise ValueError("site origin count differs from derived sites")
    groups = {source_id: [] for source_id in source_sites}
    seen = set()
    for identity, source_id, translation in result.origin_sites:
        if identity in seen or identity not in derived_sites:
            raise ValueError("site origin is incomplete or duplicated")
        seen.add(identity)
        if source_id not in source_sites:
            raise ValueError("site origin names an absent source atom")
        translation = read_integer_offset(translation)
        fractional = derived_sites[identity].fractional
        expected = matvec3(matrix, fractional)
        source_fractional = source_sites[source_id].fractional
        if any(not (0 <= x < 1) for x in fractional):
            raise ValueError("derived coordinate is not half open")
        if expected != add3(source_fractional, translation):
            raise ValueError("S y differs from f + n")
        groups[source_id].append(translation)
    if seen != set(derived_sites):
        raise ValueError("site origins do not cover derived sites")
    magnitude = abs(determinant3(matrix))
    for source_id, images in groups.items():
        if len(images) != magnitude:
            raise ValueError("source atom replica count differs from determinant")
        if len({coset_key(matrix, n) for n in images}) != magnitude:
            raise ValueError("source atom replicas repeat a coset")


def _check_ledger(source_cell, matrix, result):
    source_sites = index_sites(source_cell.sites)
    source_bonds = {bond.id: bond for bond in source_cell.bonds}
    derived_sites = index_sites(result.cell.sites)
    origin_site = {identity: (source_id, translation)
                   for identity, source_id, translation in result.origin_sites}
    origin_bond = dict(result.origin_bonds)
    if len(origin_bond) != len(result.origin_bonds):
        raise ValueError("duplicate bond origin")
    if set(origin_bond) != {row.id for row in result.bonds}:
        raise ValueError("bond origins do not cover the ledger")
    positions = dict(result.positions)
    if set(positions) != set(derived_sites):
        raise ValueError("positions do not cover derived sites")
    derived_bonds = {bond.id: bond for bond in result.cell.bonds}
    ledger_bonds = {row.id: row for row in result.bonds}
    if set(derived_bonds) != set(ledger_bonds):
        raise ValueError("derived cell bonds differ from the ledger")
    for identity, bond in derived_bonds.items():
        row = ledger_bonds[identity]
        if (row.u, row.v, row.offset) != (
                bond.u, bond.v, read_integer_offset(bond.offset)):
            raise ValueError("derived cell bond differs from the ledger")
    for identity, position in result.positions:
        expected = matvec3(result.cell.lattice.rows,
                           derived_sites[identity].fractional)
        if position != expected:
            raise ValueError("derived position differs from A S y")
    expected_volume = abs(determinant3(result.cell.lattice.rows))
    if result.volume != expected_volume:
        raise ValueError("derived volume differs from |det A S|")
    if result.cell.lattice.rows != matmul3(source_cell.lattice.rows, matrix):
        raise ValueError("derived lattice differs from A S")
    coverage = {}
    for row in result.bonds:
        source_bond = source_bonds.get(origin_bond[row.id])
        if source_bond is None:
            raise ValueError("bond origin names an absent source bond")
        u_origin = origin_site.get(row.u)
        v_origin = origin_site.get(row.v)
        if u_origin is None or v_origin is None:
            raise ValueError("bond endpoint lacks a site origin")
        if u_origin[0] != source_bond.u or v_origin[0] != source_bond.v:
            raise ValueError("bond endpoint source identity differs")
        offset = read_integer_offset(row.offset)
        key = (source_bond.id, row.u)
        if key in coverage:
            raise ValueError("source bond start replica is duplicated")
        coverage[key] = row.id
        start_y = derived_sites[row.u].fractional
        end_y = derived_sites[row.v].fractional
        delta = subtract3(add3(end_y, offset), start_y)
        expected = matvec3(source_cell.lattice.rows, subtract3(
            add3(source_sites[source_bond.v].fractional, source_bond.offset),
            source_sites[source_bond.u].fractional))
        if matvec3(result.cell.lattice.rows, delta) != expected:
            raise ValueError("derived bond displacement differs from source")
        if row.displacement != expected:
            raise ValueError("stored displacement differs from source")
        if row.length_squared != sum(
                (x * x for x in expected), Fraction(0)):
            raise ValueError("stored length squared differs from displacement")
    expected_keys = {
        (bond.id, identity)
        for bond in source_bonds.values()
        for identity, (source_id, _) in origin_site.items()
        if source_id == bond.u
    }
    if set(coverage) != expected_keys:
        raise ValueError("derived bond coverage differs from source")
    for bond_id in source_bonds:
        starts = {identity for (key_bond, identity) in coverage
                  if key_bond == bond_id}
        if len(starts) != abs(determinant3(matrix)):
            raise ValueError("source bond replica count differs from determinant")


def audit_supercell_evidence(source_cell, result):
    verify_cell_digest(source_cell, result.source[2])
    if result.source[3] != "supercell-v1":
        raise ValueError("unexpected supercell calculation")
    _check_site_origins(source_cell, result.matrix, result)
    _check_ledger(source_cell, result.matrix, result)
    return result
