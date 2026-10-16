"""SPLAS-CN 0.2: a new conservative research implementation.

Legacy Fortran is retained unmodified. This package is an independently tested
successor, not a numerical reproduction of the legacy implementation.
"""
from .model import Config, Physiology, Model, Result, simulate

__version__ = '0.2.0'
__all__ = ['Config', 'Physiology', 'Model', 'Result', 'simulate']
