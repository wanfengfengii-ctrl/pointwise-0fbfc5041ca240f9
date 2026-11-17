"""超胞派生公开契约与非答案核账夹具。

V1接口、规范cell和原回归见README.md。本题只新增如下能力。
Tasks(archive).build_supercell(product_id, revision, S) 返回整体提交的成果整数ID；
Results.read/export使用原入口。每次成功可新增成果，不要求请求去重。
S为3x3 list/tuple，行内亦为list/tuple，元素必须是非bool的Python int；
元素绝对值<=3，1<=abs(det S)<=32；源修订至多40原子、80键，空集合法。
非法参数或超出上述边界抛ValueError；合法但缺失来源抛KeyError。
沿用V1产品与修订参数规则；多个参数同时非法时不规定检查顺序。

read的超胞成果恰含source、units、matrix、cell、origins、volume、positions、bonds。
source沿用V1四字段，calculation为supercell-v1；digest仍是固定源原胞摘要。
units、volume、positions及bonds沿用V1形状、精确Fraction数值与物理定义，
但身份、晶格和偏移对应派生cell。matrix为独立的整数行列表S。
cell采用Archive.read_revision的规范nm字典格式，其中有理数仍为整数对。
origins为{sites:[{id,source_id,translation}],bonds:[{id,source_id}]}，
translation为整数三元列表；id均引用派生cell，source_id引用固定源原子或源键。
派生原子和键ID是各自非空唯一字符串，编号与枚举顺序任选；不得依赖其命名规则。

新晶格为AS。每个源原子f恰产生D=abs(det S)个副本y，0<=y_i<1，
其整数来源n满足Sy=f+n。每个源原子的n在Z^3/SZ^3中恰覆盖一次。
同坐标不同源身份不得焊接。每条源有向键在每个起点副本上恰有一条派生键，
其终点来自该源终点原子；派生位移AS(y_v+T-y_u)须精确等于源位移
A(f_v+t-f_u)。T严格为整数向量，保留每个源键ID的来源，不合并自像键、
反向键、零位移键或重复端点键。负行列式合法。不使用浮点或容差，
允许有限枚举、全量重建及复用V1几何核账和事务，不要求正规形或索引。

export返回JSON文本：format为crystal-supercell-evidence、version为1，
source与V1导出完全同形，含固定源cell和raw_cell；result与read同形，
仅将其中Fraction按V1规则编码为约分正分母整数对。
matrix、origins和派生cell均须持久保存，不能在读取或导出时重新派生。
已有V1数据库及其全部接口、结果、证据、原测试保持兼容；不指定新表或内部编码。
失败不新增本次任务、成果及来源映射，也不新增原胞修订或移动当前指针；
计算和保存异常原样传播。返回容器、S、来源后续修改以及关闭重开不得改变旧成果。
不增加空间群、猜键、最短像、并发或崩溃恢复功能。

手算例：源A=单位阵，单原子a的f=(0,3/4,0)，键b=(a,a,(0,1,0))，
S=((2,1,0),(0,1,0),(0,0,1))。两副本的(n,y)集合为
((1,0,0),(1/8,3/4,0))、((2,0,0),(5/8,3/4,0))。
从n=(1,0,0)起，终点n=(2,0,0)，T=(-1,1,0)；
从n=(2,0,0)起，终点n=(1,0,0)，T=(0,1,0)。两条位移均为(0,1,0)。
此例固定物理结果，不固定派生ID。其它合法输入仍按上述完整规则判断。

--fixtures只核对输入与手算事件，不要求V2已实现。
--verify在实现后运行全部原测试，再执行下述联合核账；--demo展示真实斜胞成果。
交付时更新.pointwise-runtime.json，使verify保留全部原回归并核验超胞，
demo实际展示斜胞发布、重开和离库核账；保留原测试字节。说明副本身份与键
偏移交接不变量，并补充独立核验；本文件是公开参照，不限制合理的额外业务用例。
"""
import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest


SKEW = [[2, 1, 0], [0, 1, 0], [0, 0, 1]]
MATRICES = [
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[2, 0, 0], [0, 1, 0], [0, 0, 1]],
    SKEW,
    [[-2, 1, 0], [0, 1, 0], [0, 0, 1]],
    [[1, 2, 0], [1, 0, 1], [0, 1, 1]],
]


def fixture():
    return {
        "unit": "nm", "lattice": [[1, [1, 2], 0], [0, 1, 0], [0, 0, 2]],
        "sites": [{"id": x, "fractional": [0, [3, 4], 0]} for x in ("a", "z")],
        "bonds": [
            {"id": "b", "u": "a", "v": "a", "offset": [0, 1, 0]},
            {"id": "dup", "u": "a", "v": "a", "offset": [0, 1, 0]},
            {"id": "back", "u": "a", "v": "a", "offset": [0, -1, 0]},
            {"id": "zero", "u": "a", "v": "z", "offset": [0, 0, 0]}],
    }


def determinant(a):
    return sum(a[0][i] * (a[1][(i+1)%3]*a[2][(i+2)%3]
                          - a[1][(i+2)%3]*a[2][(i+1)%3]) for i in range(3))


def mv(a, v):
    return [sum(a[i][j]*v[j] for j in range(3)) for i in range(3)]


def solve(a, b):
    d = determinant(a)
    return [F(determinant([[b[r] if c == k else a[r][c]
                            for c in range(3)] for r in range(3)]), d)
            for k in range(3)]


def q(pair):
    assert type(pair) is list and len(pair) == 2
    assert all(type(x) is int for x in pair) and pair[1] > 0
    v = F(*pair)
    assert [v.numerator, v.denominator] == pair
    return v


def inspect_cell(cell):
    assert cell["unit"] == "nm"
    a = [[q(x) for x in row] for row in cell["lattice"]]
    f = {s["id"]: [q(x) for x in s["fractional"]] for s in cell["sites"]}
    b = {b["id"]: b for b in cell["bonds"]}
    assert len(f) == len(cell["sites"]) and len(b) == len(cell["bonds"])
    return a, f, b


def audit(packet):
    """只消费交付证据，不枚举或生成目标超胞。"""
    assert packet["format"] == "crystal-supercell-evidence" and packet["version"] == 1
    source, r = packet["source"], packet["result"]
    canonical = json.dumps(source["cell"], ensure_ascii=True, sort_keys=True,
                           separators=(",", ":")).encode("ascii")
    assert hashlib.sha256(canonical).hexdigest() == source["digest"]
    for k in ("product_id", "revision", "digest"):
        assert r["source"][k] == source[k]
    assert r["source"]["calculation"] == "supercell-v1"
    assert r["units"] == {"length": "nm", "length_squared": "nm2", "volume": "nm3"}
    a, f, b = inspect_cell(source["cell"])
    c, y, edges = inspect_cell(r["cell"])
    s = r["matrix"]
    assert all(type(x) is int and abs(x) <= 3 for row in s for x in row)
    d = abs(determinant(s))
    assert 1 <= d <= 32
    assert c == [[sum(a[i][k]*s[k][j] for k in range(3))
                  for j in range(3)] for i in range(3)]
    assert q(r["volume"]) == abs(determinant(c))
    sites = {o["id"]: o for o in r["origins"]["sites"]}
    origins = {o["id"]: o for o in r["origins"]["bonds"]}
    assert set(sites) == set(y) and len(sites) == len(r["origins"]["sites"])
    assert set(origins) == set(edges) and len(origins) == len(r["origins"]["bonds"])
    positions = {p["id"]: p for p in r["positions"]}
    ledger = {e["id"]: e for e in r["bonds"]}
    assert set(positions) == set(y) and len(positions) == len(r["positions"])
    assert set(ledger) == set(edges) and len(ledger) == len(r["bonds"])
    assert len(y) == len(f)*d and len(edges) == len(b)*d
    groups = {identity: [] for identity in f}
    for identity, o in sites.items():
        assert o["source_id"] in f and all(0 <= x < 1 for x in y[identity])
        n = o["translation"]
        assert len(n) == 3 and all(type(x) is int for x in n)
        assert mv(s, y[identity]) == [f[o["source_id"]][i]+n[i] for i in range(3)]
        assert [q(x) for x in positions[identity]["cartesian"]] == mv(c, y[identity])
        groups[o["source_id"]].append(n)
    for images in groups.values():
        assert len(images) == d
        for i, n in enumerate(images):
            for m in images[:i]:
                assert any(x.denominator != 1 for x in solve(s, [n[j]-m[j] for j in range(3)]))
    coverage = set()
    for identity, e in edges.items():
        original = b[origins[identity]["source_id"]]
        u, v = e["u"], e["v"]
        assert sites[u]["source_id"] == original["u"]
        assert sites[v]["source_id"] == original["v"]
        key = (original["id"], u)
        assert key not in coverage
        coverage.add(key)
        row = ledger[identity]
        assert (row["u"], row["v"], row["offset"]) == (u, v, e["offset"])
        assert len(e["offset"]) == 3 and all(type(t) is int for t in e["offset"])
        delta = mv(c, [y[v][i]+e["offset"][i]-y[u][i] for i in range(3)])
        expected = mv(a, [f[original["v"]][i]+original["offset"][i]-f[original["u"]][i]
                          for i in range(3)])
        assert delta == expected == [q(x) for x in row["displacement"]]
        assert q(row["length_squared"]) == sum(x*x for x in expected)
    assert coverage == {(e["id"], identity) for e in b.values()
                        for identity, o in sites.items() if o["source_id"] == e["u"]}
    return {"sites": len(y), "bonds": len(edges), "volume": r["volume"]}


def fixtures_check():
    assert all(1 <= abs(determinant(s)) <= 32 for s in MATRICES)
    for n, y in [([1, 0, 0], [F(1, 8), F(3, 4), F(0)]),
                 ([2, 0, 0], [F(5, 8), F(3, 4), F(0)])]:
        assert mv(SKEW, y) == [n[0], F(3, 4), F(0)]
    print(json.dumps({"source": fixture(), "matrices": MATRICES,
                      "hand_case": "see module contract"}, ensure_ascii=True))


def workflow(demo=False):
    from crystal import Archive, Tasks, Results
    with tempfile.TemporaryDirectory(prefix="supercell-check-") as directory:
        path = Path(directory) / "archive.sqlite"
        retained = []
        with Archive(path) as archive:
            revision = archive.put("p", fixture())
            primitive = Tasks(archive).run_geometry("p", revision)
            primitive_packet = Results(archive).export(primitive)
            for s in ([SKEW] if demo else MATRICES):
                matrix = copy.deepcopy(s)
                rid = Tasks(archive).build_supercell("p", revision, matrix)
                matrix[0][0] = 999
                packet = json.loads(Results(archive).export(rid))
                assert packet["result"]["matrix"] == s
                assert packet["source"]["cell"] == archive.read_revision("p", revision)
                audit(packet)
                original = Results(archive).read(rid)
                value = Results(archive).read(rid)
                value["origins"]["sites"].clear()
                value["matrix"][0][0] = 999
                assert Results(archive).read(rid) == original
                retained.append((rid, original, packet))
            if not demo:
                for invalid in ([], [[1,0,0]]*3, [[4,0,0],[0,1,0],[0,0,1]],
                                [[True,0,0],[0,1,0],[0,0,1]],
                                [[F(1),0,0],[0,1,0],[0,0,1]]):
                    try:
                        Tasks(archive).build_supercell("p", revision, invalid)
                    except ValueError:
                        pass
                    else:
                        raise AssertionError("invalid matrix accepted")
                try:
                    Tasks(archive).build_supercell("absent", 1, SKEW)
                except KeyError:
                    pass
                else:
                    raise AssertionError("missing source accepted")
                before = list(archive.connection.iterdump())
                archive.connection.set_authorizer(
                    lambda action, arg1, arg2, db, context:
                    sqlite3.SQLITE_DENY if action == sqlite3.SQLITE_INSERT
                    else sqlite3.SQLITE_OK)
                try:
                    try:
                        Tasks(archive).build_supercell("p", revision, SKEW)
                    except sqlite3.DatabaseError:
                        pass
                    else:
                        raise AssertionError("database refusal did not propagate")
                finally:
                    archive.connection.set_authorizer(None)
                assert list(archive.connection.iterdump()) == before
                empty = fixture()
                empty["sites"], empty["bonds"] = [], []
                er = archive.put("empty", empty)
                audit(json.loads(Results(archive).export(
                    Tasks(archive).build_supercell("empty", er, SKEW))))
            assert archive.current_revision("p") == revision
            assert archive.connection.execute(
                "SELECT COUNT(*) FROM revisions WHERE product_id='p'").fetchone()[0] == 1
            changed = fixture()
            changed["lattice"][0][0] = 3
            archive.put("p", changed)
            assert Results(archive).export(primitive) == primitive_packet
        with Archive(path) as archive:
            for rid, original, packet in retained:
                assert Results(archive).read(rid) == original
                assert json.loads(Results(archive).export(rid)) == packet
        # 以下只用已导出的JSON，不打开数据库，不调用产品几何。
        for _, _, packet in retained:
            summary = audit(packet)
            if demo:
                print(json.dumps({"evidence": packet, "independent_audit": summary}))


def main():
    parser = argparse.ArgumentParser(description="Supercell evidence checks")
    group = parser.add_mutually_exclusive_group()
    for mode in ("fixtures", "verify", "demo"):
        group.add_argument("--" + mode, dest="mode", action="store_const", const=mode)
    mode = parser.parse_args().mode or "fixtures"
    if mode == "fixtures":
        fixtures_check()
    elif mode == "demo":
        workflow(demo=True)
    else:
        root = Path(__file__).resolve().parent
        suite = unittest.defaultTestLoader.discover(str(root/"tests"), top_level_dir=str(root))
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        if not result.wasSuccessful():
            return 1
        workflow()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
