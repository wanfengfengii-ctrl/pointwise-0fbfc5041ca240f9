"""Exact derivation of a supercell result from a fixed primitive cell.

Derivation is a pure function of the stored source revision and the integer
matrix.  Nothing here is called while reading or exporting: those paths consume
the persisted matrix, derived cell, origins and ledger instead.
"""
from .input_digest import verify_cell_digest
from .rational_matrix import matmul3
from .cell_geometry import compute_cell_geometry
from .bond_geometry import compute_bond_geometry
from .periodic_graph import build_periodic_graph
from .intake import accept_cell
from .supercell_matrix import (
    integer_fiber, read_supercell_matrix, site_representatives)
from .supercell_codec import (
    SupercellResult, decode_supercell_result, encode_supercell_result)
from .supercell_audit import audit_supercell_evidence

MAX_SOURCE_SITES = 40
MAX_SOURCE_BONDS = 80


def _check_source_size(cell):
    if len(cell.sites) > MAX_SOURCE_SITES:
        raise ValueError("source revision exceeds the site bound")
    if len(cell.bonds) > MAX_SOURCE_BONDS:
        raise ValueError("source revision exceeds the bond bound")


def _unique_identity(source_id, ordinal, used):
    """Deterministic, collision-free replica identity within one namespace."""
    candidate = f"{source_id}#{ordinal}"
    attempt = 1
    while candidate in used:
        candidate = f"{source_id}#{ordinal}#{attempt}"
        attempt += 1
    used.add(candidate)
    return candidate


def derive_supercell(cell, matrix, source):
    """Return the exact supercell result for ``cell`` and integer ``matrix``."""
    verify_cell_digest(cell, source[2])
    matrix = read_supercell_matrix(matrix)
    _check_source_size(cell)
    derived_lattice = matmul3(cell.lattice.rows, matrix)
    replicas = {}
    origin_sites = []
    derived_sites = []
    used_sites = set()
    for site in sorted(cell.sites, key=lambda item: item.id):
        images = []
        for ordinal, (translation, fractional) in enumerate(
                site_representatives(matrix, site.fractional)):
            identity = _unique_identity(site.id, ordinal, used_sites)
            images.append((translation, fractional, identity))
            derived_sites.append((identity, fractional))
            origin_sites.append((identity, site.id, translation))
        replicas[site.id] = tuple(images)
    derived_bonds = []
    origin_bonds = []
    used_bonds = set()
    for bond in sorted(cell.bonds, key=lambda item: item.id):
        for ordinal, (start_n, _, start_id) in enumerate(replicas[bond.u]):
            targets = []
            for end_n, _, end_id in replicas[bond.v]:
                translation = integer_fiber(matrix, start_n, end_n,
                                            bond.offset)
                if translation is not None:
                    targets.append((end_id, translation))
            if len(targets) != 1:
                raise ValueError("supercell bond has no unique target replica")
            end_id, translation = targets[0]
            identity = _unique_identity(bond.id, ordinal, used_bonds)
            derived_bonds.append({
                "id": identity, "u": start_id, "v": end_id,
                "offset": list(translation),
            })
            origin_bonds.append((identity, bond.id))
    derived_cell, _ = accept_cell({
        "unit": "nm",
        "lattice": [[value for value in row] for row in derived_lattice],
        "sites": [{"id": identity, "fractional": list(fractional)}
                  for identity, fractional in derived_sites],
        "bonds": derived_bonds,
    })
    grid = build_periodic_graph(derived_cell.sites, derived_cell.bonds)
    volume, positions = compute_cell_geometry(derived_cell)
    result = SupercellResult(
        source=source, matrix=matrix, cell=derived_cell,
        origin_sites=tuple(sorted(origin_sites)),
        origin_bonds=tuple(sorted(origin_bonds)),
        volume=volume, positions=positions,
        bonds=compute_bond_geometry(derived_cell.lattice, grid))
    audit_supercell_evidence(cell, result)
    return decode_supercell_result(encode_supercell_result(result))
