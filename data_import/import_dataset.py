#!/usr/bin/env python3
"""Automate descriptor-library data imports using the existing Postgres container.

The script mirrors the manual workflow documented in `data_import/directions_for_pnuema.md`
while letting each library supply its paths and column layouts through a YAML
configuration file. It runs `docker exec` commands so you do **not** need to expose
Postgres outside the container or provide connection strings.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List

try:
    import yaml
except ModuleNotFoundError as exc:  # pragma: no cover - import guard
    raise SystemExit(
        "PyYAML is required to run the data import helper. Install it with 'pip install PyYAML'."
    ) from exc


class ConfigError(RuntimeError):
    """Raised when a configuration file is invalid."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Import dataset via docker exec + psql")
    parser.add_argument(
        "--config",
        required=True,
        help="Path to YAML configuration file",
    )
    return parser.parse_args()


def load_config(path: Path) -> Dict[str, Any]:
    try:
        data = yaml.safe_load(path.read_text())
    except FileNotFoundError as exc:
        raise ConfigError(f"Configuration file not found: {path}") from exc
    except yaml.YAMLError as exc:
        raise ConfigError(f"Invalid YAML in {path}: {exc}") from exc

    if "container_name" not in data:
        raise ConfigError("Configuration missing required key 'container_name'")

    db_conf = data.get("database") or {}
    if "name" not in db_conf:
        raise ConfigError("Configuration missing database.name")

    return data


def docker_exec_psql(container: str, db: str, args: Iterable[str]) -> None:
    cmd = [
        "docker",
        "exec",
        "-i",
        container,
        "psql",
        "-v",
        "ON_ERROR_STOP=1",
        "-U",
        "postgres",
        "-d",
        db,
        *args,
    ]
    subprocess.run(cmd, check=True)


def run_sql(container: str, db: str, sql: str) -> None:
    print(f"Executing on {db}: {sql}")
    docker_exec_psql(container, db, ["-c", sql])


def run_sql_file(container: str, db: str, file_path: str) -> None:
    print(f"Applying SQL file on {db}: {file_path}")
    docker_exec_psql(container, db, ["-f", file_path])


def sql_literal(value: str) -> str:
    return value.replace("'", "''")


def quote_ident(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def qualified_identifier(name: str) -> str:
    parts = name.split(".")
    return ".".join(quote_ident(part) for part in parts)


def _table_basename(name: str) -> str:
    return name.split(".")[-1]


def resolve_table_name(table_conf: Dict[str, Any], property_sets: Dict[str, Any]) -> str:
    if "property_set" in table_conf:
        prop_name = table_conf["property_set"]
        if prop_name not in property_sets:
            raise ConfigError(f"Unknown property_set '{prop_name}'")
        prop_def = property_sets[prop_name]
        return prop_def.get("table_name", f"{prop_name}_data")
    name = table_conf.get("name")
    if not name:
        raise ConfigError("Table entry must define 'name' or 'property_set'")
    return name


def ensure_database(container: str, db_conf: Dict[str, Any]) -> None:
    db_name = db_conf["name"]
    maintenance_db = db_conf.get("maintenance_db", "postgres")
    drop_database = bool(db_conf.get("drop_database"))

    if drop_database:
        print(f"Dropping database {db_name} if it exists")
        terminate_sql = (
            "SELECT pg_terminate_backend(pid) FROM pg_stat_activity "
            f"WHERE datname = '{sql_literal(db_name)}' AND pid <> pg_backend_pid();"
        )
        run_sql(container, maintenance_db, terminate_sql)
        run_sql(container, maintenance_db, f"DROP DATABASE IF EXISTS {quote_ident(db_name)};")

    create_sql = (
        "DO $$BEGIN IF NOT EXISTS (SELECT FROM pg_database WHERE datname = '"
        + sql_literal(db_name)
        + "') THEN CREATE DATABASE "
        + quote_ident(db_name)
        + "; END IF; END$$;"
    )
    run_sql(container, maintenance_db, create_sql)

    owner = db_conf.get("owner")
    if owner:
        run_sql(
            container,
            maintenance_db,
            f"ALTER DATABASE {quote_ident(db_name)} OWNER TO {quote_ident(owner)};",
        )


def apply_extensions(container: str, db_name: str, extensions: Iterable[str]) -> None:
    for ext in extensions:
        run_sql(
            container,
            db_name,
            f"CREATE EXTENSION IF NOT EXISTS {quote_ident(ext)};",
        )


def run_statements(container: str, db_name: str, statements: Iterable[str], stage: str) -> None:
    statements = [stmt.strip() for stmt in statements or [] if stmt and stmt.strip()]
    if not statements:
        return
    print(f"Running {stage} statements ({len(statements)} commands)")
    for stmt in statements:
        run_sql(container, db_name, stmt)


def ensure_property_table(
    container: str,
    db_name: str,
    table_name: str,
    definition: Dict[str, Any],
) -> List[str]:
    columns = definition.get("columns")
    if not columns:
        raise ConfigError(f"Property set for table '{table_name}' must define columns")

    create_lines: List[str] = []
    copy_columns: List[str] = []

    for col in columns:
        if isinstance(col, str):
            column_name = col
            column_type = "text"
            nullable = True
            default = None
            csv_name = column_name
        elif isinstance(col, dict):
            if "name" not in col or "type" not in col:
                raise ConfigError(
                    f"Property column definitions require 'name' and 'type' (table {table_name})."
                )
            column_name = col["name"]
            column_type = col["type"]
            nullable = col.get("nullable", True)
            default = col.get("default")
            csv_name = col.get("csv", column_name)
        else:
            raise ConfigError(
                "Property column definitions must be strings or mapping objects with name/type."
            )

        column_sql = f"{quote_ident(column_name)} {column_type}"
        if not nullable:
            column_sql += " NOT NULL"
        if default is not None:
            column_sql += f" DEFAULT {default}"
        create_lines.append(column_sql)
        copy_columns.append(column_name)

    primary_key = definition.get("primary_key") or []
    if primary_key:
        pk_cols = ", ".join(quote_ident(col) for col in primary_key)
        create_lines.append(f"PRIMARY KEY ({pk_cols})")

    create_sql = (
        f"CREATE TABLE IF NOT EXISTS {qualified_identifier(table_name)} (\n    "
        + ",\n    ".join(create_lines)
        + "\n);"
    )
    run_sql(container, db_name, create_sql)

    # Optional molecule_id foreign key
    molecule_fk_enabled = definition.get(
        "molecule_fk",
        any(col_name == "molecule_id" for col_name in copy_columns),
    )
    if molecule_fk_enabled and "molecule_id" in copy_columns:
        constraint_name = f"{_table_basename(table_name)}_molecule_id_fkey"
        table_regclass = f"'{table_name}'::regclass"
        fk_sql = (
            "DO $$BEGIN IF NOT EXISTS ("
            "SELECT 1 FROM pg_constraint WHERE conname = '{conname}' AND conrelid = {regclass}) "
            "THEN ALTER TABLE {table} ADD CONSTRAINT {con_ident} FOREIGN KEY (molecule_id) "
            "REFERENCES public.molecule(molecule_id){on_delete}; END IF; END$$;"
        ).format(
            conname=sql_literal(constraint_name),
            regclass=table_regclass,
            table=qualified_identifier(table_name),
            con_ident=quote_ident(constraint_name),
            on_delete=f" ON DELETE {definition.get('molecule_fk_on_delete')}" if definition.get("molecule_fk_on_delete") else "",
        )
        run_sql(container, db_name, fk_sql)

    for raw_index_sql in definition.get("indexes", []) or []:
        index_sql = raw_index_sql.format(table=qualified_identifier(table_name))
        run_sql(container, db_name, index_sql)

    for extra_sql in definition.get("post_create_sql", []) or []:
        run_sql(container, db_name, extra_sql.format(table=qualified_identifier(table_name)))

    return copy_columns


def _normalize_columns(columns: List[Any], table_name: str) -> List[str]:
    result: List[str] = []
    for col in columns:
        if isinstance(col, str):
            result.append(col)
        elif isinstance(col, dict) and "name" in col:
            result.append(col["name"])
        else:
            raise ConfigError(
                f"Invalid column entry in table '{table_name}'. Use string names or objects with 'name'."
            )
    return result


def copy_into_table(container: str, db_name: str, table_conf: Dict[str, Any]) -> None:
    table_name = table_conf.get("name")
    property_set_name = table_conf.get("property_set")

    if property_set_name:
        property_sets = table_conf["__property_sets__"]
        if property_set_name not in property_sets:
            raise ConfigError(f"Unknown property_set '{property_set_name}'")
        prop_def = property_sets[property_set_name]
        table_name = prop_def.get("table_name", f"{property_set_name}_data")
        copy_columns = ensure_property_table(
            container,
            db_name,
            table_name,
            prop_def,
        )
    else:
        if not table_name:
            raise ConfigError("Each table entry must include a 'name' or 'property_set'")
        columns = table_conf.get("columns") or []
        if not columns:
            raise ConfigError(f"Table '{table_name}' must define a non-empty 'columns' list")
        copy_columns = _normalize_columns(columns, table_name)

    csv_path = table_conf.get("csv_path")
    if not csv_path:
        raise ConfigError(f"Table '{table_name}' must define 'csv_path'")

    copy_options = table_conf.get("copy_options", "FORMAT csv, HEADER true")

    column_list = ", ".join(quote_ident(col) for col in copy_columns)
    escaped_path = csv_path.replace("'", "''")
    copy_sql = (
        f"\\copy {qualified_identifier(table_name)} ({column_list}) FROM '{escaped_path}' WITH ({copy_options});"
    )

    print(f"Copying data into {table_name} from {csv_path}")
    docker_exec_psql(container, db_name, ["-c", copy_sql])


def process_tables(
    container: str,
    db_name: str,
    tables: List[Dict[str, Any]],
    property_sets: Dict[str, Any],
) -> None:
    for table in tables or []:
        pre_sql = table.get("pre_sql")
        post_sql = table.get("post_sql")

        target_table_name = resolve_table_name(table, property_sets)

        # Inject property_set definitions when requested
        table_for_copy = table
        if "property_set" in table:
            table_for_copy = {**table, "__property_sets__": property_sets}

        if pre_sql:
            formatted_pre = [stmt.format(table=target_table_name) for stmt in pre_sql]
            run_statements(container, db_name, formatted_pre, f"pre-import for {target_table_name}")

        copy_into_table(container, db_name, table_for_copy)

        if post_sql:
            formatted_sql = [stmt.format(table=target_table_name) for stmt in post_sql]
            run_statements(container, db_name, formatted_sql, f"post-import for {target_table_name}")


def main() -> None:
    args = parse_args()
    config_path = Path(args.config).resolve()

    try:
        config = load_config(config_path)
    except ConfigError as exc:
        print(f"Configuration error: {exc}", file=sys.stderr)
        sys.exit(1)

    container = config["container_name"]
    db_conf = config["database"]
    db_name = db_conf["name"]

    try:
        ensure_database(container, db_conf)
        apply_extensions(container, db_name, config.get("extensions", []))

        run_statements(container, db_name, config.get("pre_schema_sql"), "pre-schema")

        for schema_file in config.get("schema_files", []):
            run_sql_file(container, db_name, schema_file)

        run_statements(container, db_name, config.get("post_schema_sql"), "post-schema")

        process_tables(
            container,
            db_name,
            config.get("tables", []),
            config.get("property_sets", {}) or {},
        )

        run_statements(container, db_name, config.get("post_import_sql"), "post-import")
    except subprocess.CalledProcessError as exc:
        print(f"Command failed with exit code {exc.returncode}", file=sys.stderr)
        sys.exit(exc.returncode or 2)
    except ConfigError as exc:
        print(f"Configuration error: {exc}", file=sys.stderr)
        sys.exit(1)

    print("Import completed successfully.")


if __name__ == "__main__":
    main()
