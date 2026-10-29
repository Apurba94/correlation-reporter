# correlation-reporter

**Publication-ready correlation tables: Pearson or Spearman matrices with significance stars,
p-values, pairwise n and Fisher-z confidence intervals, from any CSV.**

Part of the *Research Automation Toolkit*: 100 open tools for daily research work.

## Why

"Table 1: Means, SDs and intercorrelations" appears in countless papers. Building it by hand from
statistics software output is tedious and error-prone (APA formatting, lower triangle, stars,
missing data). This tool produces it directly, using pairwise-complete observations.

## Features

- Pearson or Spearman (average ranks for ties), pairwise-complete n
- Exact two-sided p-values from the t distribution (incomplete beta), Fisher-z 95% CIs
- APA-style Markdown matrix: numbered variables, lower triangle, no leading zeros, true minus, stars
- Long CSV output (a, b, r, p, n, CI) for further processing or multiple-testing correction
- Non-numeric columns are skipped automatically

## Usage

```bash
pip install .
correlation-reporter survey.csv -c age,income,wellbeing,stress
correlation-reporter data.csv -m spearman -f long > correlations.csv
```

## Testing

```bash
python -m unittest discover -s tests -t .
```

## License

MIT (c) JANIN A APURBA
