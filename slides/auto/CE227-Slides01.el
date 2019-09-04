(TeX-add-style-hook
 "CE227-Slides01"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("babel" "brazil") ("inputenc" "utf8") ("caption" "hang") ("lato" "default") ("abntex2cite" "alf" "abnt-emphasize=bf" "abnt-etal-list=2" "abnt-and-type=&")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-environments-local "semiverbatim")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "header"
    "beamer"
    "beamer10"
    "babel"
    "inputenc"
    "tikz"
    "graphicx"
    "upquote"
    "listings"
    "hyperref"
    "amsmath"
    "color"
    "caption"
    "lato"
    "inconsolata"
    "abntex2cite")
   (TeX-add-symbols
    '("tcDR" 1)
    '("tcDB" 1)
    '("blue" 1)
    '("red" 1))
   (LaTeX-add-bibliographies
    "references.bib")
   (LaTeX-add-xcolor-definecolors
    "DarkBlue"))
 :latex)

