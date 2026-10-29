"""Correlation matrices with p-values, ready for papers."""
from .core import pearson, spearman, corr_pvalue, correlation_matrix, to_markdown

__version__ = "0.1.0"
__all__ = ["pearson", "spearman", "corr_pvalue", "correlation_matrix", "to_markdown"]
