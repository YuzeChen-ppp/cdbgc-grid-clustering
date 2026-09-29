# LaTeX source of the poster

`poster.tex` is the A0 landscape research poster for this project
(*Cell Density-Based Grid Clustering*), typeset with the **Imperial College London
poster template**.

## Not included in this repository

The `ImperialPoster` class file and the Imperial fonts/logos are **not** redistributed here.
The template is © Imperial College London and is licensed for use by Imperial staff and
students for university business only. Get it from the official source:

- https://github.com/ImperialCollegeLondon/imperial_latex_templates

Drop `ImperialPoster.cls` (and the template's supporting files) next to `poster.tex` before
compiling.

## Compiling

The template requires **XeLaTeX**, not pdfLaTeX, because of the custom fonts:

```bash
xelatex poster.tex
bibtex poster
xelatex poster.tex
xelatex poster.tex
```

(the three passes are for the table of contents-free float placement and the bibliography).

Expected structure when compiling:

```
latex/
├── poster.tex
├── ImperialPoster.cls          <- bring your own (official template)
├── siam_1.9781611973440.96.bib
└── Images/
    ├── Grid construction.png
    ├── heat map.png
    ├── grid.png
    ├── DBCV.png
    └── comparison.png
```

## Notes

- `siam_1.9781611973440.96.bib` was **reconstructed** from the reference printed on the
  poster (Moulavi et al., *Density-Based Clustering Validation*, SDM 2014). Swap in your
  original file if it differs.
- A compiled copy of the poster is at `../poster/CDBGC_poster.pdf`.
