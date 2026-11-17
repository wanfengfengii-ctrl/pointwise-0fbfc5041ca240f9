import sys
sys.dont_write_bytecode = True
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
from crystal import Archive, Tasks, Results

SUPERCELL = [[2, 1, 0], [0, 1, 0], [0, 0, 1]]


def determinant(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


def offline_values(packet):
    cell = packet["source"]["cell"]
    a = [[F(*x) for x in row] for row in cell["lattice"]]
    f = {s["id"]:[F(*x) for x in s["fractional"]] for s in cell["sites"]}
    rows = {}
    for bond in cell["bonds"]:
        delta = [f[bond["v"]][j]+bond["offset"][j]-f[bond["u"]][j]
                 for j in range(3)]
        d = [sum(a[i][j]*delta[j] for j in range(3)) for i in range(3)]
        rows[bond["id"]] = {"displacement":d,
                           "length_squared":sum(x*x for x in d)}
    return abs(determinant(a)),rows


def offline_supercell(packet):
    """Recompute a supercell packet from exported JSON only.

    The identity handover is checked against the exported ``origins``: each
    source key must appear on exactly |det S| start replicas, each derived
    displacement must equal the source displacement, and every replica's
    integer translation must reproduce ``S y = f + n``.
    """
    source, result = packet["source"], packet["result"]
    a = [[F(*x) for x in row] for row in source["cell"]["lattice"]]
    f = {s["id"]:[F(*x) for x in s["fractional"]]
         for s in source["cell"]["sites"]}
    b = {e["id"]:e for e in source["cell"]["bonds"]}
    c = [[F(*x) for x in row] for row in result["cell"]["lattice"]]
    y = {s["id"]:[F(*x) for x in s["fractional"]]
         for s in result["cell"]["sites"]}
    s = result["matrix"]
    d = abs(determinant(s))
    assert c == [[sum(a[i][k]*s[k][j] for k in range(3)) for j in range(3)]
                 for i in range(3)]
    origins = {o["id"]:o for o in result["origins"]["sites"]}
    bond_origins = {o["id"]:o for o in result["origins"]["bonds"]}
    starts = {}
    for identity, origin in origins.items():
        n = origin["translation"]
        assert [sum(s[i][j]*y[identity][j] for j in range(3))
                for i in range(3)] == [
            f[origin["source_id"]][i]+n[i] for i in range(3)]
        starts.setdefault((origin["source_id"], tuple(n)), []).append(identity)
    assert all(len(v) == 1 for v in starts.values())
    per_bond = {}
    for row in result["bonds"]:
        original = b[bond_origins[row["id"]]["source_id"]]
        assert origins[row["u"]]["source_id"] == original["u"]
        assert origins[row["v"]]["source_id"] == original["v"]
        per_bond[original["id"]] = per_bond.get(original["id"], 0) + 1
        delta = [y[row["v"]][j]+row["offset"][j]-y[row["u"]][j]
                 for j in range(3)]
        expected = [sum(a[i][j]*(f[original["v"]][j]+original["offset"][j]
                                 -f[original["u"]][j]) for j in range(3))
                    for i in range(3)]
        assert [sum(c[i][j]*delta[j] for j in range(3))
                for i in range(3)] == expected
        assert [F(*x) for x in row["displacement"]] == expected
        assert F(*row["length_squared"]) == sum(x*x for x in expected)
    assert set(per_bond.values()) == {d}
    return {"det": d, "sites": len(y), "bonds": len(result["bonds"])}


def demonstrate(path, product, revision, sample=None):
    with Archive(path) as archive:
        if sample is not None:
            revision = archive.put(product,sample)
        rid = Tasks(archive).run_geometry(product,revision)
        old = Results(archive).read(rid)
        if sample is not None:
            changed = json.loads(json.dumps(sample))
            changed["lattice"][0][0] = 2
            new_revision = archive.put(product,changed)
            newer = Results(archive).read(
                Tasks(archive).run_geometry(product,new_revision))
            assert newer["volume"] == F(4)
            assert old["volume"] == F(2)
    with Archive(path) as archive:
        assert Results(archive).read(rid) == old
        packet = json.loads(Results(archive).export(rid))
    volume, expected = offline_values(packet)
    assert old["volume"] == volume
    for row in old["bonds"]:
        assert row["displacement"] == expected[row["id"]]["displacement"]
        assert row["length_squared"] == expected[row["id"]]["length_squared"]
    if sample is not None:
        assert expected["uv"]["displacement"] == [F(5,4),F(3,2),F(0)]
        assert expected["uv"]["length_squared"] == F(61,16)
    def pair(q):
        return [q.numerator,q.denominator]
    expected_json = {
        key:{"displacement":[pair(x) for x in row["displacement"]],
             "length_squared":pair(row["length_squared"])}
        for key,row in expected.items()}
    actual_json = {
        row["id"]:{"displacement":row["displacement"],
                    "length_squared":row["length_squared"]}
        for row in packet["result"]["bonds"]}
    scenario = {
        "case_id":"N8",
        "input":packet["source"]["cell"],
        "expected":{"volume":pair(volume),"bonds":expected_json},
        "actual":{"volume":packet["result"]["volume"],"bonds":actual_json},
        "downstream_expected":{k:v["length_squared"]
                               for k,v in expected_json.items()},
        "downstream_actual":{k:v["length_squared"]
                             for k,v in actual_json.items()},
        "code_paths":["crystal/archive.py","crystal/tasks.py",
                      "crystal/executor.py","crystal/result_store.py",
                      "crystal/results.py","crystal/export.py"],
        "oracle":"Fraction recomputation from exported source after archive close"}
    print("POINTWISE_SCENARIO " + json.dumps(scenario))
    return 0


def demonstrate_supercell(path, product, revision):
    """Publish a skew supercell, close, reopen and audit outside the archive."""
    with Archive(path) as archive:
        rid = Tasks(archive).build_supercell(product, revision, SUPERCELL)
        original = Results(archive).read(rid)
    with Archive(path) as archive:
        assert Results(archive).read(rid) == original
        packet = json.loads(Results(archive).export(rid))
    summary = offline_supercell(packet)
    print("POINTWISE_SCENARIO " + json.dumps({
        "case_id":"N9",
        "input":{"matrix":SUPERCELL,"source":packet["source"]["cell"]},
        "expected":{"det":2,"sites":2*summary["det"],
                    "bonds":len(packet["source"]["cell"]["bonds"])
                    *summary["det"]},
        "actual":{"volume":packet["result"]["volume"],
                  "origins_sites":len(packet["result"]["origins"]["sites"]),
                  "origins_bonds":len(packet["result"]["origins"]["bonds"]),
                  "independent_audit":summary},
        "downstream_expected":[4,1],
        "downstream_actual":packet["result"]["volume"],
        "code_paths":["crystal/tasks.py","crystal/supercell.py",
                      "crystal/supercell_matrix.py","crystal/supercell_codec.py",
                      "crystal/supercell_audit.py","crystal/result_store.py",
                      "crystal/results.py","crystal/export.py"],
        "oracle":"exact coset and displacement replay from exported JSON "
                 "after archive close"}))
    return 0


def main():
    parser = argparse.ArgumentParser(description="Primitive-cell laboratory workflow")
    parser.add_argument("--database",type=Path)
    parser.add_argument("--product")
    parser.add_argument("--revision",type=int)
    args = parser.parse_args()
    if args.database is not None:
        if not args.database.is_file() or not args.product or args.revision is None:
            parser.error("existing --database requires --product and --revision")
        return demonstrate_supercell(args.database,args.product,args.revision)
    if args.product is not None or args.revision is not None:
        parser.error("--product and --revision require --database")
    sample = json.loads((Path(__file__).resolve().parent/
                        "examples/day_cell.json").read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="crystal-demo-") as directory:
        path = Path(directory)/"lab.sqlite"
        status = demonstrate(path,"day-cell",None,sample)
        if status:
            return status
        return demonstrate_supercell(path,"day-cell",1)


if __name__ == "__main__":
    raise SystemExit(main())
