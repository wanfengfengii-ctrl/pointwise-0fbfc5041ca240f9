from dataclasses import dataclass
from .result_codec import CALCULATION, SUPERCELL_CALCULATION
from .executor import execute_geometry
from .supercell import derive_supercell
from .supercell_matrix import read_supercell_matrix
from .result_store import ResultStore


@dataclass(frozen=True)
class GeometryTask:
    product_id: str
    revision: int
    digest: str
    calculation: str = CALCULATION

    @property
    def source(self):
        return (self.product_id, self.revision, self.digest, self.calculation)


@dataclass(frozen=True)
class SupercellTask:
    product_id: str
    revision: int
    digest: str
    calculation: str = SUPERCELL_CALCULATION

    @property
    def source(self):
        return (self.product_id, self.revision, self.digest, self.calculation)


class Tasks:
    def __init__(self, archive):
        self.revisions = archive.revisions
        self.store = ResultStore(archive.connection)

    def run_geometry(self, product_id, revision):
        fixed = self.revisions.read_revision(product_id, revision)
        task = GeometryTask(fixed.product_id, fixed.revision, fixed.digest)
        result = execute_geometry(fixed.cell, task.source)
        return self.store.publish_bundle(task, result)

    def build_supercell(self, product_id, revision, S):
        matrix = read_supercell_matrix(S)
        fixed = self.revisions.read_revision(product_id, revision)
        task = SupercellTask(fixed.product_id, fixed.revision, fixed.digest)
        result = derive_supercell(fixed.cell, matrix, task.source)
        return self.store.publish_supercell_bundle(task, result)
