(TeX-add-style-hook
 "ce227-av01"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt" "a4")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("babel" "brazil" "brazilian") ("geometry" "top=1cm" "left=1cm" "bottom=1cm" "right=1cm" "nohead")))
   (TeX-run-style-hooks
    "latex2e"
    "formulas1-prob"
    "article"
    "art11"
    "inputenc"
    "babel"
    "geometry"
    "amsmath"
    "graphicx"
    "comment"
    "array"
    "Sweave")
   (LaTeX-add-labels
    "fig:vero-e-vero-aprox"
    "fig:gamma-inv"
    "fig:sim-gamainv"))
 :latex)

