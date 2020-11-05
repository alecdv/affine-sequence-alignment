# Affine global sequence alignment

The program takes two sample DNA sequences and outputs the highest scoring
global alignment using the specified scoring scheme and affine gap penalty.

Input file contains three lines: scoring values, two DNA strings to align.
Scoring values are:

1. p1 - value for match
2. p2 - penlaty for mismatch
3. g - penlaty for opening gap
4. s - penalty for extending gap

Example:
    5 3 1 2
    BCDA
    ABCDEFG

Scoring scheme for alignment of two bases x, y:
    S(x, y) = { p1 if x = y, -p2 if x != y }

Affine gap penalty for gap of length k:
    w(k) = { -g - s(k - 1) if k >= 1, 0 if k = 0}
