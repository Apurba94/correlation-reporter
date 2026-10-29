import io
import os
import tempfile
import unittest
from contextlib import redirect_stdout

from correlation_reporter.cli import main
from correlation_reporter.core import correlation_matrix, corr_pvalue, fisher_ci, pearson, spearman, to_markdown


class CorrTests(unittest.TestCase):
    def test_pearson_and_spearman(self):
        x = [1, 2, 3, 4, 5]
        self.assertAlmostEqual(pearson(x, [2, 4, 6, 8, 10]), 1.0)
        self.assertAlmostEqual(pearson(x, [5, 4, 3, 2, 1]), -1.0)
        self.assertAlmostEqual(spearman(x, [1, 4, 9, 16, 100]), 1.0)     # monotone but not linear
        self.assertLess(pearson(x, [1, 4, 9, 16, 100]), 1.0)

    def test_pvalue_reference(self):
        # r = .5, n = 20 -> t(18) = 2.449, p = .0247 (tabulated)
        self.assertAlmostEqual(corr_pvalue(0.5, 20), 0.0247, places=3)
        self.assertAlmostEqual(corr_pvalue(0.0, 30), 1.0)
        self.assertEqual(corr_pvalue(1.0, 10), 0.0)

    def test_fisher_ci(self):
        lo, hi = fisher_ci(0.5, 20)
        self.assertAlmostEqual(lo, 0.069, places=2)
        self.assertAlmostEqual(hi, 0.77, places=2)

    def test_matrix_pairwise_complete(self):
        names, res = correlation_matrix({"a": [1, 2, 3, 4, None], "b": [2, 4, 6, 8, 10], "c": [1, None, 2, None, 3]})
        self.assertEqual(res[("a", "b")][2], 4)
        self.assertEqual(res[("a", "c")][2], 2)
        self.assertNotEqual(res[("a", "c")][0], res[("a", "c")][0])   # NaN when n < 3
        md = to_markdown(names, res)
        self.assertIn("| (2) b | 1.00*** | — |  |", md)

    def test_cli(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "d.csv")
            with open(p, "w", encoding="utf-8") as fh:
                fh.write("x,y,label\n1,2,a\n2,1,b\n3,4,c\n4,3,d\n5,6,e\n")
            with redirect_stdout(io.StringIO()) as buf:
                main([p, "-f", "long", "-m", "spearman"])
            lines = buf.getvalue().splitlines()
            self.assertEqual(len(lines), 2)
            self.assertTrue(lines[1].startswith("x,y,0.8000"))


if __name__ == "__main__":
    unittest.main()
