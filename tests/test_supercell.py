from fractions import Fraction as F
import copy
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from crystal import Archive, Tasks, Results
from crystal.rational_matrix import determinant3
from crystal.supercell import derive_supercell
from crystal.supercell_codec import encode_supercell_result, decode_supercell_result
from crystal.supercell_matrix import read_supercell_matrix

SKEW = [[2, 1, 0], [0, 1, 0], [0, 0, 1]]
MATRICES = [
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[2, 0, 0], [0, 1, 0], [0, 0, 1]],
    SKEW,
    [[-2, 1, 0], [0, 1, 0], [0, 0, 1]],
    [[1, 2, 0], [1, 0, 1], [0, 1, 1]],
]


def hand_cell():
    return {"unit": "nm", "lattice": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            "sites": [{"id": "a", "fractional": [0, F(3, 4), 0]}],
            "bonds": [{"id": "b", "u": "a", "v": "a",
                       "offset": [0, 1, 0]}]}


def rich_cell():
    return {
        "unit": "nm", "lattice": [[1, F(1, 2), 0], [0, 1, 0], [0, 0, 2]],
        "sites": [{"id": "a", "fractional": [0, F(3, 4), 0]},
                  {"id": "z", "fractional": [0, F(3, 4), 0]}],
        "bonds": [
            {"id": "b", "u": "a", "v": "a", "offset": [0, 1, 0]},
            {"id": "dup", "u": "a", "v": "a", "offset": [0, 1, 0]},
            {"id": "back", "u": "a", "v": "a", "offset": [0, -1, 0]},
            {"id": "zero", "u": "a", "v": "z", "offset": [0, 0, 0]}],
    }


def q(pair):
    assert type(pair) is list and len(pair) == 2
    assert all(type(x) is int for x in pair) and pair[1] > 0
    value = F(*pair)
    assert [value.numerator, value.denominator] == pair
    return value


def mv(matrix, vector):
    return [sum(matrix[i][j] * vector[j] for j in range(3)) for i in range(3)]


def replay(packet):
    """Independent replay that only consumes exported JSON.

    It checks the lattice relation, the half-open replica coverage by cosets
    and the displacement handover; it never calls the product geometry.
    """
    def determinant(a):
        return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
                -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
                +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))

    assert packet["format"] == "crystal-supercell-evidence"
    assert packet["version"] == 1
    source, result = packet["source"], packet["result"]
    canonical = json.dumps(source["cell"], sort_keys=True,
                           separators=(",", ":"), ensure_ascii=True)
    assert hashlib.sha256(canonical.encode("ascii")).hexdigest() == source["digest"]
    for key in ("product_id", "revision", "digest"):
        assert result["source"][key] == source[key]
    assert result["source"]["calculation"] == "supercell-v1"
    a = [[q(x) for x in row] for row in source["cell"]["lattice"]]
    f = {s["id"]: [q(x) for x in s["fractional"]]
         for s in source["cell"]["sites"]}
    b = {e["id"]: e for e in source["cell"]["bonds"]}
    c = [[q(x) for x in row] for row in result["cell"]["lattice"]]
    y = {s["id"]: [q(x) for x in s["fractional"]]
         for s in result["cell"]["sites"]}
    edges = {e["id"]: e for e in result["cell"]["bonds"]}
    s = result["matrix"]
    d = abs(determinant(s))
    assert 1 <= d <= 32
    assert c == [[sum(a[i][k]*s[k][j] for k in range(3)) for j in range(3)]
                 for i in range(3)]
    assert q(result["volume"]) == abs(determinant(c))
    origins = {o["id"]: o for o in result["origins"]["sites"]}
    bond_origins = {o["id"]: o for o in result["origins"]["bonds"]}
    assert set(origins) == set(y) == {p["id"] for p in result["positions"]}
    assert set(bond_origins) == set(edges) == {r["id"] for r in result["bonds"]}
    positions = {p["id"]: p for p in result["positions"]}
    ledger = {r["id"]: r for r in result["bonds"]}
    groups = {identity: [] for identity in f}
    for identity, origin in origins.items():
        assert origin["source_id"] in f
        assert all(0 <= x < 1 for x in y[identity])
        n = origin["translation"]
        assert len(n) == 3 and all(type(x) is int for x in n)
        assert mv(s, y[identity]) == [
            f[origin["source_id"]][i] + n[i] for i in range(3)]
        assert [q(x) for x in positions[identity]["cartesian"]] == mv(c, y[identity])
        groups[origin["source_id"]].append(n)
    for images in groups.values():
        assert len(images) == d

    def adjugate(matrix):
        def minor(i, j):
            rows = [r for r in range(3) if r != i]
            cols = [c for c in range(3) if c != j]
            return (matrix[rows[0]][cols[0]]*matrix[rows[1]][cols[1]]
                    - matrix[rows[0]][cols[1]]*matrix[rows[1]][cols[0]])
        cof = [[(-1)**(i+j)*minor(i, j) for j in range(3)] for i in range(3)]
        return [[cof[j][i] for j in range(3)] for i in range(3)]

    sign = 1 if determinant(s) > 0 else -1
    adj = adjugate(s)

    def class_key(n):
        return tuple((sign*sum(adj[j][k]*n[k] for k in range(3))) % d
                     for j in range(3))

    for images in groups.values():
        assert len({class_key(n) for n in images}) == d
    coverage = set()
    for identity, edge in edges.items():
        original = b[bond_origins[identity]["source_id"]]
        u, v = edge["u"], edge["v"]
        assert origins[u]["source_id"] == original["u"]
        assert origins[v]["source_id"] == original["v"]
        assert (original["id"], u) not in coverage
        coverage.add((original["id"], u))
        offset = edge["offset"]
        assert len(offset) == 3 and all(type(x) is int for x in offset)
        row = ledger[identity]
        assert (row["u"], row["v"], row["offset"]) == (u, v, offset)
        delta = mv(c, [y[v][i]+offset[i]-y[u][i] for i in range(3)])
        expected = mv(a, [f[original["v"]][i]+original["offset"][i]
                          - f[original["u"]][i] for i in range(3)])
        assert delta == expected == [q(x) for x in row["displacement"]]
        assert q(row["length_squared"]) == sum(
            (x*x for x in expected), F(0))
    assert coverage == {(e["id"], identity) for e in b.values()
                        for identity, origin in origins.items()
                        if origin["source_id"] == e["u"]}
    return len(y), len(edges)


class TestSupercellHandCase(unittest.TestCase):
    def test_hand_example_exact(self):
        with Archive() as archive:
            revision = archive.put("p", hand_cell())
            rid = Tasks(archive).build_supercell("p", revision, SKEW)
            value = Results(archive).read(rid)
        self.assertEqual(value["matrix"], SKEW)
        self.assertEqual(value["volume"], F(2))
        sites = {o["id"]: o for o in value["origins"]["sites"]}
        self.assertEqual(len(sites), 2)
        self.assertEqual({o["translation"][0] for o in sites.values()},
                         {1, 2})
        fractions = {s["id"]: [q(x) for x in s["fractional"]]
                     for s in value["cell"]["sites"]}
        for identity, origin in sites.items():
            self.assertEqual(origin["source_id"], "a")
            self.assertTrue(all(0 <= x < 1 for x in fractions[identity]))
            n = origin["translation"]
            self.assertEqual(
                [2*fractions[identity][0] + fractions[identity][1],
                 fractions[identity][1], fractions[identity][2]],
                [n[0], F(3, 4) + n[1], n[2]])
        rows = {b["id"]: b for b in value["bonds"]}
        self.assertEqual(len(rows), 2)
        self.assertEqual({tuple(b["offset"]) for b in rows.values()},
                         {(-1, 1, 0), (0, 1, 0)})
        for row in rows.values():
            self.assertEqual(row["displacement"], [F(0), F(1), F(0)])
            self.assertEqual(row["length_squared"], F(1))
        # Each start replica gets exactly one bond and distinct replica ends.
        starts = {row["u"] for row in rows.values()}
        ends = {row["v"] for row in rows.values()}
        self.assertEqual(starts, ends)
        self.assertEqual(starts, set(sites))

    def test_skewed_lattice_handover(self):
        with Archive() as archive:
            revision = archive.put("p", rich_cell())
            primitive = Tasks(archive).run_geometry("p", revision)
            for matrix in MATRICES:
                given = copy.deepcopy(matrix)
                rid = Tasks(archive).build_supercell("p", revision, given)
                given[0][0] = 999
                packet = json.loads(Results(archive).export(rid))
                self.assertEqual(packet["result"]["matrix"], matrix)
                self.assertEqual(packet["source"]["cell"],
                                 archive.read_revision("p", revision))
                sites, bonds = replay(packet)
                self.assertEqual(sites, 2 * abs(determinant3(matrix)))
                self.assertEqual(bonds, 4 * abs(determinant3(matrix)))
            self.assertEqual(Results(archive).read(primitive)["source"]
                             ["calculation"], "primitive-v1")


class TestSupercellPersistence(unittest.TestCase):
    def test_frozen_across_close_and_correction(self):
        with tempfile.TemporaryDirectory(prefix="supercell-test-") as directory:
            path = Path(directory) / "archive.sqlite"
            with Archive(path) as archive:
                revision = archive.put("p", rich_cell())
                primitive = Tasks(archive).run_geometry("p", revision)
                primitive_packet = Results(archive).export(primitive)
                rid = Tasks(archive).build_supercell("p", revision, SKEW)
                original = Results(archive).read(rid)
                packet = json.loads(Results(archive).export(rid))
                changed = rich_cell()
                changed["lattice"][0][0] = 3
                newer = archive.put("p", changed)
                Tasks(archive).build_supercell("p", newer, SKEW)
                self.assertEqual(archive.current_revision("p"), newer)
                self.assertEqual(archive.connection.execute(
                    "SELECT COUNT(*) FROM revisions WHERE product_id='p'"
                ).fetchone()[0], 2)
                self.assertEqual(Results(archive).export(primitive),
                                 primitive_packet)
                mutable = Results(archive).read(rid)
                mutable["matrix"][0][0] = 999
                mutable["origins"]["sites"].clear()
                mutable["bonds"].clear()
                self.assertEqual(Results(archive).read(rid), original)
            with Archive(path) as archive:
                self.assertEqual(Results(archive).read(rid), original)
                self.assertEqual(json.loads(Results(archive).export(rid)), packet)
                self.assertEqual(Results(archive).export(primitive),
                                 primitive_packet)
                replay(packet)
                self.assertEqual(archive.connection.execute(
                    "PRAGMA foreign_key_check").fetchall(), [])

    def test_repeated_publication_does_not_require_dedup(self):
        with Archive() as archive:
            revision = archive.put("p", rich_cell())
            first = Tasks(archive).build_supercell("p", revision, SKEW)
            second = Tasks(archive).build_supercell("p", revision,
                                                    [list(row) for row in SKEW])
            self.assertNotEqual(first, second)
            self.assertEqual(Results(archive).read(first),
                             Results(archive).read(second))

    def test_failure_rolls_back_and_keeps_sources(self):
        with Archive() as archive:
            revision = archive.put("p", rich_cell())
            before_rows = tuple(archive.connection.execute(
                "SELECT COUNT(*) FROM " + table).fetchone()[0]
                for table in ("tasks", "results", "result_bodies", "revisions"))
            before = list(archive.connection.iterdump())
            archive.connection.execute(
                "CREATE TEMP TRIGGER reject_body BEFORE INSERT ON result_bodies "
                "BEGIN SELECT RAISE(ABORT,'body write rejected'); END")
            with self.assertRaisesRegex(sqlite3.IntegrityError,
                                        "body write rejected"):
                Tasks(archive).build_supercell("p", revision, SKEW)
            archive.connection.execute("DROP TRIGGER reject_body")
            self.assertEqual(list(archive.connection.iterdump()), before)
            after_rows = tuple(archive.connection.execute(
                "SELECT COUNT(*) FROM " + table).fetchone()[0]
                for table in ("tasks", "results", "result_bodies", "revisions"))
            self.assertEqual(after_rows, before_rows)
            self.assertEqual(archive.current_revision("p"), revision)

    def test_invalid_matrices_and_missing_source(self):
        with Archive() as archive:
            revision = archive.put("p", rich_cell())
            before = tuple(archive.connection.execute(
                "SELECT COUNT(*) FROM " + table).fetchone()[0]
                for table in ("tasks", "results", "result_bodies"))
            invalid = [
                [], [[1, 0, 0]] * 2, [[1, 0, 0]] * 4,
                [[4, 0, 0], [0, 1, 0], [0, 0, 1]],
                [[True, 0, 0], [0, 1, 0], [0, 0, 1]],
                [[F(1), 0, 0], [0, 1, 0], [0, 0, 1]],
                [[1.0, 0, 0], [0, 1, 0], [0, 0, 1]],
                [[0, 0, 0], [0, 1, 0], [0, 0, 1]],
                [[2, 0, 0], [0, 4, 0], [0, 0, 4]],
                [[1, 0], [0, 1, 0], [0, 0, 1]],
                "identity",
            ]
            for matrix in invalid:
                with self.subTest(matrix=matrix), self.assertRaises(ValueError):
                    Tasks(archive).build_supercell("p", revision, matrix)
            with self.assertRaises(KeyError):
                Tasks(archive).build_supercell("absent", 1, SKEW)
            with self.assertRaises(KeyError):
                Tasks(archive).build_supercell("p", revision + 1, SKEW)
            with self.assertRaises(ValueError):
                Tasks(archive).build_supercell("p", 0, SKEW)
            self.assertEqual(tuple(archive.connection.execute(
                "SELECT COUNT(*) FROM " + table).fetchone()[0]
                for table in ("tasks", "results", "result_bodies")), before)

    def test_empty_source(self):
        with Archive() as archive:
            cell = rich_cell()
            cell["sites"], cell["bonds"] = [], []
            revision = archive.put("empty", cell)
            rid = Tasks(archive).build_supercell("empty", revision, SKEW)
            packet = json.loads(Results(archive).export(rid))
        self.assertEqual(packet["result"]["positions"], [])
        self.assertEqual(packet["result"]["bonds"], [])
        self.assertEqual(packet["result"]["origins"],
                         {"sites": [], "bonds": []})
        self.assertEqual(replay(packet), (0, 0))


class TestSupercellBoundsAndMigration(unittest.TestCase):
    def test_source_size_bounds(self):
        with Archive() as archive:
            sites = [{"id": f"s{i}", "fractional": [F(i, 41), 0, 0]}
                     for i in range(41)]
            bonds = [{"id": f"b{i}", "u": "s0", "v": "s0",
                      "offset": [i, 0, 0]} for i in range(81)]
            over = {"unit": "nm", "lattice": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    "sites": sites, "bonds": bonds[:1]}
            revision = archive.put("over-sites", over)
            with self.assertRaises(ValueError):
                Tasks(archive).build_supercell("over-sites", revision, SKEW)
            over["sites"], over["bonds"] = sites[:40], bonds
            revision = archive.put("over-bonds", over)
            with self.assertRaises(ValueError):
                Tasks(archive).build_supercell("over-bonds", revision, SKEW)
            over["bonds"] = bonds[:80]
            revision = archive.put("at-bound", over)
            rid = Tasks(archive).build_supercell("at-bound", revision, SKEW)
            value = Results(archive).read(rid)
            self.assertEqual(len(value["positions"]), 80)
            self.assertEqual(len(value["bonds"]), 160)
            self.assertEqual(value["source"]["calculation"], "supercell-v1")

    def test_large_determinant_handover(self):
        matrix = [[-3, 3, 3], [3, -3, 3], [3, 3, -3]]
        self.assertEqual(determinant3(matrix), 108)
        with Archive() as archive:
            revision = archive.put("p", rich_cell())
            with self.assertRaises(ValueError):
                Tasks(archive).build_supercell("p", revision, matrix)
            matrix = [[-3, 1, 0], [0, -2, 1], [0, 0, -3]]
            self.assertEqual(determinant3(matrix), -18)
            rid = Tasks(archive).build_supercell("p", revision, matrix)
            packet = json.loads(Results(archive).export(rid))
        self.assertEqual(replay(packet), (2 * 18, 4 * 18))

    def test_v1_database_gains_supercell_support(self):
        with tempfile.TemporaryDirectory(prefix="supercell-v1-") as directory:
            path = Path(directory) / "archive.sqlite"
            with Archive(path) as archive:
                revision = archive.put("p", rich_cell())
                primitive = Tasks(archive).run_geometry("p", revision)
                old_value = Results(archive).read(primitive)
                old_packet = Results(archive).export(primitive)
            connection = sqlite3.connect(path, isolation_level=None)
            connection.execute("PRAGMA foreign_keys=OFF")
            connection.execute("BEGIN IMMEDIATE")
            connection.execute(
                "CREATE TABLE old_tasks (task_id INTEGER PRIMARY KEY "
                "AUTOINCREMENT, product_id TEXT NOT NULL, revision INTEGER "
                "NOT NULL, digest TEXT NOT NULL, calculation TEXT NOT NULL "
                "CHECK(calculation = 'primitive-v1'), FOREIGN KEY("
                "product_id,revision) REFERENCES revisions(product_id,revision))")
            connection.execute("INSERT INTO old_tasks SELECT * FROM tasks")
            connection.execute("DROP TABLE tasks")
            connection.execute("ALTER TABLE old_tasks RENAME TO tasks")
            connection.execute("COMMIT")
            connection.close()
            with Archive(path) as archive:
                self.assertEqual(Results(archive).read(primitive), old_value)
                self.assertEqual(Results(archive).export(primitive), old_packet)
                rid = Tasks(archive).build_supercell("p", revision, SKEW)
                self.assertGreater(rid, primitive)
                replay(json.loads(Results(archive).export(rid)))
                self.assertEqual(archive.connection.execute(
                    "PRAGMA foreign_key_check").fetchall(), [])
                self.assertEqual(archive.current_revision("p"), revision)


class TestSupercellCodecAndAudit(unittest.TestCase):
    def test_codec_roundtrip_and_negative_determinant(self):
        cell, digest = self._accepted(rich_cell())
        for matrix in MATRICES:
            with self.subTest(matrix=matrix):
                result = derive_supercell(
                    cell, matrix, ("p", 1, digest, "supercell-v1"))
                encoded = encode_supercell_result(result)
                self.assertEqual(decode_supercell_result(encoded), result)
                value = json.loads(encoded)
                del value["result"]["origins"]["sites"]
                with self.assertRaises(ValueError):
                    decode_supercell_result(json.dumps(value))
                value = json.loads(encoded)
                value["result"]["origins"]["sites"].pop()
                with self.assertRaises(ValueError):
                    decode_supercell_result(json.dumps(value))
                value = json.loads(encoded)
                value["result"]["matrix"][0][0] = True
                with self.assertRaises(ValueError):
                    decode_supercell_result(json.dumps(value))

    def test_stored_body_tampering_is_rejected(self):
        with Archive() as archive:
            revision = archive.put("p", rich_cell())
            rid = Tasks(archive).build_supercell("p", revision, SKEW)
            payload = archive.connection.execute(
                "SELECT payload FROM result_bodies WHERE result_id=?",
                (rid,)).fetchone()[0]
            tampered = json.loads(payload)
            tampered["result"]["bonds"][0]["length_squared"] = [99, 1]
            archive.connection.execute(
                "UPDATE result_bodies SET payload=? WHERE result_id=?",
                (json.dumps(tampered), rid))
            with self.assertRaises(ValueError):
                Results(archive).read(rid)
            with self.assertRaises(ValueError):
                Results(archive).export(rid)

    def test_matrix_and_nullspace_semantics(self):
        # A singular or oversized matrix is refused before any derivation.
        for matrix in ([[0, 0, 0], [0, 1, 0], [0, 0, 1]],
                       [[2, 0, 0], [0, 2, 0], [0, 0, 9]],
                       [[3, 3, 3], [3, 3, 3], [3, 3, 3]]):
            with self.subTest(matrix=matrix), self.assertRaises(ValueError):
                read_supercell_matrix(matrix)
        self.assertEqual(determinant3(read_supercell_matrix(
            [[-2, 1, 0], [0, 1, 0], [0, 0, 1]])), -2)

    def _accepted(self, value):
        from crystal.intake import accept_cell
        return accept_cell(value)


if __name__ == "__main__":
    unittest.main()
