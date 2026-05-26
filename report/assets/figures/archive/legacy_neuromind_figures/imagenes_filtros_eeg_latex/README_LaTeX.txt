# Imágenes generadas para LaTeX: filtros EEG

Carpeta creada a partir de los conceptos de las infografías proporcionadas, sin recortar ni reutilizar elementos de las imágenes originales.

Formatos incluidos por cada figura:
- PDF vectorial: recomendado para LaTeX.
- SVG vectorial: útil para edición posterior en Inkscape/Illustrator.
- PNG a 300 dpi: útil para previsualización o inserción directa.

Sugerencia LaTeX:
\begin{figure}[htbp]
  \centering
  \includegraphics[width=0.85\textwidth]{01_paso_alto_ordenes.pdf}
  \caption{Respuesta en frecuencia de un filtro paso alto Butterworth para distintos órdenes.}
\end{figure}

Listado:
01_paso_alto_ordenes
02_paso_bajo_ordenes
03_notch_50hz_respuesta
04_bandas_eeg
05_resumen_tres_filtros
06_fase_cero_filtfilt
07_pipeline_procesamiento_eeg
08_paso_alto_antes_despues
09_paso_bajo_antes_despues
10_notch_antes_despues
11_fir_vs_iir
12_regla_muestreo
