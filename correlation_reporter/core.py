"""Pairwise-complete Pearson/Spearman correlations, t-based p-values and Fisher-z confidence intervals."""
import math
import statistics


def _betacf(a, b, x):
    tiny, qab, qap, qam = 1e-300, a + b, a + 1, a - 1
    c, d = 1.0, 1 - qab * x / qap
    d = 1 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        for aa in (m * (b - m) * x / ((qam + m2) * (a + m2)), -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))):
            d = 1 + aa * d
            d = 1 / (d if abs(d) > tiny else tiny)
            c = 1 + aa / c
            c = c if abs(c) > tiny else tiny
            h *= d * c
        if abs(d * c - 1) < 3e-14:
            break
    return h


def _betainc(a, b, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbt = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x)
    if x < (a + 1) / (a + b + 2):
        return math.exp(lbt) * _betacf(a, b, x) / a
    return 1 - math.exp(lbt) * _betacf(b, a, 1 - x) / b


def corr_pvalue(r, n):
    """Two-sided p-value for H0: rho = 0 using t = r sqrt((n-2)/(1-r^2)) with n-2 df."""
    if n < 3:
        return float("nan")
    if abs(r) >= 1:
        return 0.0
    df = n - 2
    t = r * math.sqrt(df / (1 - r * r))
    return _betainc(df / 2, 0.5, df / (df + t * t))


def fisher_ci(r, n, level=0.95):
    if n < 4 or abs(r) >= 1:
        return (r, r)
    z = math.atanh(r)
    se = 1 / math.sqrt(n - 3)
    q = statistics.NormalDist().inv_cdf(1 - (1 - level) / 2)
    return (math.tanh(z - q * se), math.tanh(z + q * se))


def pearson(x, y):
    mx, my = statistics.fmean(x), statistics.fmean(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sx = math.sqrt(sum((a - mx) ** 2 for a in x))
    sy = math.sqrt(sum((b - my) ** 2 for b in y))
    return sxy / (sx * sy) if sx and sy else float("nan")


def _rank(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    ranks = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return ranks


def spearman(x, y):
    return pearson(_rank(x), _rank(y))


def correlation_matrix(columns, method="pearson"):
    """columns: dict name -> list of floats/None. Pairwise-complete. Returns dict (a, b) -> (r, p, n, ci)."""
    f = pearson if method == "pearson" else spearman
    names = list(columns)
    res = {}
    for i, a in enumerate(names):
        for b in names[i:]:
            pairs = [(x, y) for x, y in zip(columns[a], columns[b]) if x is not None and y is not None]
            n = len(pairs)
            if a == b:
                res[(a, b)] = (1.0, 0.0, n, (1.0, 1.0))
                continue
            if n < 3:
                res[(a, b)] = res[(b, a)] = (float("nan"), float("nan"), n, (float("nan"),) * 2)
                continue
            r = f([p[0] for p in pairs], [p[1] for p in pairs])
            out = (r, corr_pvalue(r, n), n, fisher_ci(r, n))
            res[(a, b)] = res[(b, a)] = out
    return names, res


def stars(p):
    if p != p:
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def to_markdown(names, res, lower=True, decimals=2):
    head = "| | " + " | ".join("(%d)" % (i + 1) for i in range(len(names))) + " |"
    lines = [head, "|" + "---|" * (len(names) + 1)]
    for i, a in enumerate(names):
        cells = []
        for j, b in enumerate(names):
            if lower and j > i:
                cells.append("")
            elif i == j:
                cells.append("—")
            else:
                r, p, n, _ = res[(a, b)]
                cells.append("" if r != r else ("%.*f" % (decimals, r)).replace("0.", ".", 1).replace("-.", "−.") + stars(p))
        lines.append("| (%d) %s | %s |" % (i + 1, a, " | ".join(cells)))
    lines.append("")
    lines.append("*Note.* * p < .05, ** p < .01, *** p < .001 (two-sided); pairwise-complete n.")
    return "\n".join(lines) + "\n"
