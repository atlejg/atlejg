"""
alerts exported via:
   gh api repos/equinor/apo-r/dependabot/alerts > github.alerts
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from pprint import pprint
from typing import Any


def load_alerts(path: str | Path) -> list[dict[str, Any]]:
    with open(path, encoding="utf-16") as f:
        return json.load(f)


def package_name(alert: dict[str, Any]) -> str:
    try:
        return alert["dependency"]["package"]["name"]
    except KeyError:
        return "<unknown>"


def dependency_scope(alert: dict[str, Any]) -> str:
    dependency = alert.get("dependency", {})
    return dependency.get("scope", "<unknown>")


def dependency_path(alert: dict[str, Any]) -> list:
    """
    Try several commonly-used locations where GitHub might expose
    dependency path information.

    This is intentionally defensive because the schema changes over time.
    """
    candidates = [
        ("dependency", "dependency_path"),
        ("dependency", "path"),
        ("dependency_path",),
        ("path",),
    ]

    for candidate in candidates:
        node: Any = alert

        try:
            for key in candidate:
                node = node[key]

            if isinstance(node, list):
                return [str(x) for x in node]

        except (KeyError, TypeError):
            pass

    return []


def root_dependency(alert: dict[str, Any]) -> str | None:
    """
    Given:

        apo-r -> fluidmagic-toolkit -> urllib3

    return:

        fluidmagic-toolkit
    """
    path = dependency_path(alert)

    if len(path) >= 2:
        return path[1]

    return None


def main() -> None:
    alerts = load_alerts("github.alerts")

    print(f"Loaded {len(alerts)} alerts\n")

    package_counter: Counter[str] = Counter()
    scope_counter: Counter[str] = Counter()
    root_counter: Counter[str] = Counter()

    for alert in alerts:
        package_counter[package_name(alert)] += 1
        scope_counter[dependency_scope(alert)] += 1

        root = root_dependency(alert)
        if root is not None:
            root_counter[root] += 1

    print("=" * 80)
    print("Top vulnerable packages")
    print("=" * 80)

    for package, count in package_counter.most_common():
        print(f"{count:4d}  {package}")

    print()

    print("=" * 80)
    print("Dependency scopes")
    print("=" * 80)

    for scope, count in scope_counter.most_common():
        print(f"{count:4d}  {scope}")

    print()

    if root_counter:
        print("=" * 80)
        print("Root dependencies introducing alerts")
        print("=" * 80)

        for root, count in root_counter.most_common():
            print(f"{count:4d}  {root}")

        print()
    else:
        print("=" * 80)
        print("No dependency-path information found")
        print("=" * 80)
        print()

    print("=" * 80)
    print("First alert structure")
    print("=" * 80)

    pprint(alerts[0])

if __name__ == '__main__':

    main()
