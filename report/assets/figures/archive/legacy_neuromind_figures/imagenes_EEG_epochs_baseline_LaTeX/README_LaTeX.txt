# Imágenes EEG para LaTeX

Carpeta generada con imágenes nuevas e independientes inspiradas en los conceptos de las infografías aportadas, sin recortar ni reutilizar partes de las imágenes originales.

Contenido:
- PNG a 300 dpi, adecuados para inclusión directa.
- PDF vectorial, recomendado para LaTeX por máxima nitidez.
- SVG vectorial, útil para edición posterior.

Sugerencia en LaTeX:
\begin{figure}[htbp]
  \centering
  \includegraphics[width=0.85\textwidth]{01_senal_EEG_continua_resting_state.pdf}
  \caption{Señal EEG continua en resting-state.}
\end{figure}
