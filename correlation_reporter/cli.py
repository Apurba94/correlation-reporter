"""Command-line interface for correlation-reporter."""
import argparse
import csv
import sys

from .core import correlation_matrix, to_markdown


def main(argv=None):
    ap = argparse.ArgumentParser(prog="correlation-reporter", description="Correlation matrix with p-values.")
    ap.add_argument("csvfile")
    ap.add_argument("-c", "--columns", help="comma-separated columns (default: all numeric)")
    ap.add_argument("-m", "--method", choices=["pearson", "spearman"], default="pearson")
    ap.add_argument("-f", "--format", choices=["markdown", "long"], default="markdown",
                    help="markdown matrix or long CSV (a, b, r, p, n, ci)")
    ap.add_argument("--full", action="store_true", help="show the full matrix instead of the lower triangle")
    args = ap.parse_args(argv)
    with open(args.csvfile, encoding="utf-8-sig", newline="") as fh:
        r = csv.DictReader(fh)
        rows, fields = list(r), r.fieldnames
    wanted = [c.strip() for c in args.columns.split(",")] if args.columns else fields
    cols = {}
    for c in wanted:
        vals = []
        for row in rows:
            try:
                vals.append(float(row[c]))
            except (ValueError, TypeError):
                vals.append(None)
        if sum(v is not None for v in vals) >= 3:
            cols[c] = vals
    names, res = correlation_matrix(cols, args.method)
    if args.format == "long":
        w = csv.writer(sys.stdout)
        w.writerow(["a", "b", "r", "p", "n", "ci_low", "ci_high"])
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                r_, p, n, ci = res[(a, b)]
                w.writerow([a, b, "%.4f" % r_, "%.4g" % p, n, "%.4f" % ci[0], "%.4f" % ci[1]])
    else:
        sys.stdout.write(to_markdown(names, res, lower=not args.full))
    return 0
