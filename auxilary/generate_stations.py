#!/usr/bin/env python3
"""Generate a WaveQLab3D station list block.

Writes a file containing:

    !---begin:station_list---
    x y z
    ...
    !---end:station_list---

Coordinates are written in Fortran-style double precision literals (e.g. 12.0d0).

Usage:
    generate_stations.py x=[start,end,step] y=[start,end,step] z=[start,end,step] OUTPUT_FILE

Examples:
    # Single station at the origin
    ./generate_stations.py x=0 y=0 z=0 stations.txt

    # Line of stations along x from 0 to 10 every 1
    ./generate_stations.py x=[0,10,1] y=0 z=0 stations.txt

    # 2D grid (x-y plane) at z=5
    ./generate_stations.py x=[0,2,1] y=[0,3,1] z=5 stations.txt

Notes:
    - Ranges are inclusive (end is included when it lands on a step).
    - If step is 0, start must equal end.
"""

import argparse
import sys
from pathlib import Path


USAGE = (
        "Usage: generate_stations.py x=[start,end,step] y=[start,end,step] z=[start,end,step] OUTPUT_FILE\n"
        "Example: ./generate_stations.py x=[0,10,1] y=0 z=0 stations.txt\n"
)


def parse_range(spec: str):
    spec = spec.strip()
    if spec.startswith("[") and spec.endswith("]"):
        inner = spec[1:-1]
        parts = [p.strip() for p in inner.split(",") if p.strip()]
        if len(parts) != 3:
            raise ValueError(f"range specification must have 3 numbers: {spec}")
        start, end, step = map(float, parts)
        if abs(step) < 1e-12:
            if abs(start - end) > 1e-12:
                raise ValueError("step cannot be zero unless start == end")
            return [start]
        values = []
        current = start
        limit = end + (1e-12 if step > 0 else -1e-12)
        while (step > 0 and current <= limit) or (step < 0 and current >= limit):
            values.append(current)
            current += step
        return values
    return [float(spec)]


def format_value(value: float) -> str:
    if abs(value - round(value)) < 1e-12:
        return f"{int(round(value))}.0d0"
    text = f"{value:.10f}".rstrip("0").rstrip(".")
    return f"{text}d0"


def load_dimension(arg_list, key, default="[0,0,0]"):
    for token in arg_list:
        if token.startswith(f"{key}="):
            return token.split("=", 1)[1]
    return default


def main():
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        print(__doc__.rstrip())
        sys.exit(0)

    if len(sys.argv) < 2:
        print(USAGE, file=sys.stderr)
        print("Run with --help to see more examples.", file=sys.stderr)
        sys.exit(1)

    args = sys.argv[1:]
    output_path = args[-1]
    specs = args[:-1]

    x_spec = load_dimension(specs, "x")
    y_spec = load_dimension(specs, "y")
    z_spec = load_dimension(specs, "z")

    x_values = parse_range(x_spec)
    y_values = parse_range(y_spec)
    z_values = parse_range(z_spec)

    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    with out_file.open("w", encoding="utf-8") as fh:
        fh.write("!---begin:station_list---\n")
        for x in x_values:
            for y in y_values:
                for z in z_values:
                    fh.write(
                        f"{format_value(x)} {format_value(y)} {format_value(z)}\n"
                    )
        fh.write("!---end:station_list---\n")


if __name__ == "__main__":
    main()
