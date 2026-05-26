# Configuracion principal de latexmk (pdfLaTeX + Biber + SyncTeX).
# Auxiliares en build/pdf.
$pdf_mode = 1;
$out_dir = 'build/pdf';
$pdflatex = 'pdflatex -interaction=nonstopmode -synctex=1 %O %S';
$biber = 'biber %O %B';
$max_repeat = 8;
