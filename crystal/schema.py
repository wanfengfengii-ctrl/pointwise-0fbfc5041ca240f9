SCHEMA_VERSION = 1
CALCULATION_CHECK = "calculation IN ('primitive-v1', 'supercell-v1')"

TABLES = (
    """CREATE TABLE IF NOT EXISTS products (
        product_id TEXT PRIMARY KEY,
        current_revision INTEGER,
        FOREIGN KEY(product_id,current_revision)
            REFERENCES revisions(product_id,revision)
            DEFERRABLE INITIALLY DEFERRED)""",
    """CREATE TABLE IF NOT EXISTS revisions (
        product_id TEXT NOT NULL REFERENCES products(product_id),
        revision INTEGER NOT NULL CHECK(revision > 0),
        digest TEXT NOT NULL,
        canonical_json TEXT NOT NULL,
        raw_json TEXT NOT NULL,
        PRIMARY KEY(product_id,revision))""",
    """CREATE TABLE IF NOT EXISTS tasks (
        task_id INTEGER PRIMARY KEY AUTOINCREMENT,
        product_id TEXT NOT NULL,
        revision INTEGER NOT NULL,
        digest TEXT NOT NULL,
        calculation TEXT NOT NULL CHECK(""" + CALCULATION_CHECK + """),
        FOREIGN KEY(product_id,revision)
            REFERENCES revisions(product_id,revision))""",
    """CREATE TABLE IF NOT EXISTS results (
        result_id INTEGER PRIMARY KEY AUTOINCREMENT,
        task_id INTEGER NOT NULL UNIQUE REFERENCES tasks(task_id))""",
    """CREATE TABLE IF NOT EXISTS result_bodies (
        result_id INTEGER PRIMARY KEY REFERENCES results(result_id),
        payload TEXT NOT NULL)""",
)


def check_schema_version(connection):
    version = connection.execute("PRAGMA user_version").fetchone()[0]
    if version not in (0, SCHEMA_VERSION):
        raise ValueError("incompatible database schema")
    if version == 0:
        tables = {row[0] for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table' "
            "AND name NOT LIKE 'sqlite_%'")}
        if tables:
            raise ValueError("unversioned nonempty database")
    return version


def _tasks_table_sql(connection):
    row = connection.execute(
        "SELECT sql FROM sqlite_master WHERE type='table' AND name='tasks'"
    ).fetchone()
    return None if row is None else row[0]


def _migrate_tasks_calculation(connection):
    """Widen the persisted calculation CHECK without changing any row.

    V1 databases only accept primitive-v1 in ``tasks``.  The table is rebuilt
    with the same columns, same identifiers and the same foreign key, while
    foreign key enforcement is suspended (SQLite would otherwise refuse to drop
    a referenced table).  The rebuild is atomic and keeps task ids, the
    AUTOINCREMENT sequence, result rows and bodies identical.
    """
    connection.execute("PRAGMA foreign_keys=OFF")
    try:
        connection.execute("BEGIN IMMEDIATE")
        try:
            connection.execute(
                """CREATE TABLE tasks_new (
                    task_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    product_id TEXT NOT NULL,
                    revision INTEGER NOT NULL,
                    digest TEXT NOT NULL,
                    calculation TEXT NOT NULL CHECK("""
                + CALCULATION_CHECK + """),
                    FOREIGN KEY(product_id,revision)
                        REFERENCES revisions(product_id,revision))""")
            connection.execute(
                "INSERT INTO tasks_new(task_id,product_id,revision,digest,"
                "calculation) SELECT task_id,product_id,revision,digest,"
                "calculation FROM tasks")
            connection.execute("DROP TABLE tasks")
            connection.execute("ALTER TABLE tasks_new RENAME TO tasks")
            connection.execute("COMMIT")
        except BaseException:
            if connection.in_transaction:
                connection.execute("ROLLBACK")
            raise
    finally:
        connection.execute("PRAGMA foreign_keys=ON")


def _needs_calculation_migration(connection):
    sql = _tasks_table_sql(connection)
    return sql is not None and "supercell-v1" not in sql


def initialize_schema(connection):
    check_schema_version(connection)
    if _needs_calculation_migration(connection):
        _migrate_tasks_calculation(connection)
    connection.execute("BEGIN IMMEDIATE")
    try:
        for statement in TABLES:
            connection.execute(statement)
        connection.execute("PRAGMA user_version=1")
        connection.execute("COMMIT")
    except BaseException:
        if connection.in_transaction:
            connection.execute("ROLLBACK")
        raise
