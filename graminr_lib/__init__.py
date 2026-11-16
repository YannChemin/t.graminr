"""
MODULE:     t.graminr (graminr_lib)

AUTHOR(S):  Yann Chemin

PURPOSE:    Grassland models of the GraMinR package, ported to Python and
            vectorised over the pixels of a GRASS GIS region.

            GraMinR (Grassland Models in R) gathers several grassland
            models behind one interface; this package ports five of them
            -- LINGRA, ModVege, BASGRA, GrasProg and GrasSim. PaSim is not
            ported: it wraps an external executable that GraMinR itself
            does not distribute.

SPDX-License-Identifier: Unlicense
"""

__all__ = [
    "basgra",
    "core",
    "cutdates",
    "fetch",
    "gio",
    "grasprog",
    "grassim",
    "lingra",
    "modvege",
    "params",
]
