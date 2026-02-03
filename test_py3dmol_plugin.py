#!/usr/bin/env python3
"""
Test script to check if genice3-py3dmol plugin is visible to GenIce3.
"""

import sys
import os
from importlib.metadata import entry_points

print("=" * 60)
print("Testing genice3-py3dmol plugin visibility")
print("=" * 60)
print(f"\nPython version: {sys.version}")
print(f"Python executable: {sys.executable}")
print(f"Python path:")
for p in sys.path:
    print(f"  {p}")

# Check entry points
print("\n1. Checking entry points for 'genice3_exporter':")
try:
    eps = list(entry_points(group="genice3_exporter"))
    print(f"   Total entry points found: {len(eps)}")
    found = False
    for ep in eps:
        print(
            f"   Found: {ep.name} -> {ep.value} (from {ep.dist.name if hasattr(ep, 'dist') else 'unknown'})"
        )
        if ep.name == "py3dmol":
            found = True
            print(f"   ✓ py3dmol entry point found!")
            try:
                module = ep.load()
                print(f"   ✓ Module loaded successfully: {module}")
                print(
                    f"   ✓ Module location: {module.__file__ if hasattr(module, '__file__') else 'N/A'}"
                )
            except Exception as e:
                print(f"   ✗ Failed to load module: {e}")
                import traceback

                traceback.print_exc()
    if not found:
        print("   ✗ py3dmol entry point NOT found")
        print(
            "   Hint: Make sure genice3-py3dmol is installed with: pip install -e /path/to/genice3-py3dmol"
        )
except Exception as e:
    print(f"   ✗ Error checking entry points: {e}")
    import traceback

    traceback.print_exc()

# Check if module can be imported directly
print("\n2. Checking direct import:")
try:
    import genice3_py3dmol

    print(f"   ✓ genice3_py3dmol imported: {genice3_py3dmol}")
    print(f"   ✓ Version: {genice3_py3dmol.__version__}")
except ImportError as e:
    print(f"   ✗ Failed to import genice3_py3dmol: {e}")

# Check if exporter module can be imported
print("\n3. Checking exporter module import:")
try:
    from genice3_py3dmol.exporter import py3dmol

    print(f"   ✓ genice3_py3dmol.exporter.py3dmol imported: {py3dmol}")
    print(f"   ✓ Has 'get_view' function: {hasattr(py3dmol, 'get_view')}")
    print(f"   ✓ Has 'desc' dict: {hasattr(py3dmol, 'desc')}")
except ImportError as e:
    print(f"   ✗ Failed to import exporter module: {e}")

# Check GenIce3 plugin system
print("\n4. Checking GenIce3 plugin system:")
try:
    from genice3.plugin import safe_import

    print("   ✓ genice3.plugin.safe_import imported")
    try:
        module = safe_import("exporter", "py3dmol")
        print(f"   ✓ safe_import('exporter', 'py3dmol') succeeded: {module}")
        print(f"   ✓ Module has 'get_view': {hasattr(module, 'get_view')}")
    except Exception as e:
        print(f"   ✗ safe_import('exporter', 'py3dmol') failed: {e}")
        import traceback

        traceback.print_exc()
except ImportError as e:
    print(f"   ✗ Failed to import genice3.plugin: {e}")

# Check Exporter function
print("\n5. Checking Exporter function:")
try:
    from genice3.plugin import Exporter

    print("   ✓ genice3.plugin.Exporter imported")
    try:
        exporter_module = Exporter("py3dmol")
        print(f"   ✓ Exporter('py3dmol') succeeded: {exporter_module}")
        print(f"   ✓ Module has 'get_view': {hasattr(exporter_module, 'get_view')}")
    except Exception as e:
        print(f"   ✗ Exporter('py3dmol') failed: {e}")
        import traceback

        traceback.print_exc()
except ImportError as e:
    print(f"   ✗ Failed to import Exporter: {e}")

print("\n" + "=" * 60)






