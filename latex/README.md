# LaTeX source of the poster

`poster.tex` is the A0 landscape research poster for this project
(*Cell Density-Based Grid Clustering*), typeset with the **Imperial College London
research poster template** (XeLaTeX).

> In the original Overleaf project the main file is named `template_poster_landscape.tex`
> (the poster was written directly into the template file rather than renamed).
> It is committed here as `poster.tex` for clarity; the body is unchanged.

## What is deliberately NOT in this repository

The template itself is © Imperial College London, 2024 — "for use of Imperial staff and
students only for university business. All rights reserved to the copyright owners."
Its assets are therefore **not redistributed** here:

| Excluded file (present in the Overleaf export) | Why |
|---|---|
| `ImperialPoster.cls` | template class file |
| `ImperialLetterhead.cls` | template class file |
| `Fonts/ImperialSansText-*.ttf` (7 weights) | Imperial's licensed typeface |
| `Images/ICL_Logo_Blue_RGB.pdf`, `Images/ICL_Logo_Blue_CMYK.pdf` | Imperial logo |
| `Images/UKRI_logo.png` | third-party (UKRI) logo |
| `template_poster_portrait.tex`, `README.md`, `LICENSE.txt` | template distribution files |

Get the official template from:

- https://github.com/ImperialCollegeLondon/imperial_latex_templates

## Compiling

The template requires **XeLaTeX**, not pdfLaTeX, because of the custom fonts:

```bash
xelatex poster.tex
bibtex poster
xelatex poster.tex
xelatex poster.tex
```

Expected layout when compiling (files marked *bring your own* come from the official template):

```
latex/
├── poster.tex
├── ImperialPoster.cls                        * bring your own
├── Fonts/ImperialSansText-*.ttf              * bring your own
├── Images/ICL_Logo_Blue_RGB.pdf              * bring your own
├── siam_1.9781611973440.96.bib
└── Images/                                   (all figures below are in this repo)
    ├── Grid construction.png                 referenced by the poster
    ├── heat map.png                          referenced by the poster
    ├── grid.png                              referenced by the poster
    ├── DBCV.png                              referenced by the poster
    ├── comparison.png                        referenced by the poster
    ├── Kmeans.png                            in the Overleaf project, not placed on the poster
    ├── CDBGC image.png                       in the Overleaf project, not placed on the poster
    └── workflow.png                          in the Overleaf project, not placed on the poster
```

## Notes

- `siam_1.9781611973440.96.bib` is the original bibliography file from the Overleaf project.
- A compiled copy of the poster is at `../poster/CDBGC_poster.pdf`.
