FUNCTION fetch_layer, layer, step, label

  ip1             = scan_parameters('ip1', 0, label)
  ip2             = scan_parameters('ip2', 0, label)
   n3             = scan_parameters( 'n3', 0, label)

  n_lines         = LONG64(0)
  
  datafile        = fetch_datafile(layer, step, label)

  IF (STRCMP(datafile, 'file_not_found')) THEN BEGIN
    PRINT, "fetch_layer: ERROR - cannot open file for layer ", layer, " of step ", step
    A             = DBLARR(1,1)
  ENDIF ELSE BEGIN
;   PRINT, "fetch_layer: opening file ", datafile, "..."
    IF (STRMATCH(datafile, '*.gz', /FOLD_CASE) EQ 1) THEN BEGIN
      OPENR, data_unit, datafile, /COMPRESS, /GET_LUN, ERROR = op_err
      n_lines     = FILE_LINES(datafile, /COMPRESS)
    ENDIF ELSE BEGIN
      OPENR, data_unit, datafile, /GET_LUN, ERROR = op_err
      n_lines     = FILE_LINES(datafile)
    ENDELSE
  ENDELSE

  local_layer     =  ((layer - 1) mod n3 ) + 1

; PRINT, "fetch_layer: layer       = ", layer
; PRINT, "fetch_layer: local_layer = ", local_layer

  lines_per_layer = n_lines / n3

  l_start         = (local_layer - 1) * lines_per_layer + 1
  l_stop          = l_start + lines_per_layer - 1
 
  IF (op_err NE 0) THEN BEGIN
    PRINT, "fetch_layer: ERROR - file exists but cannot be opened"
    A             = DBLARR(1,1)
  ENDIF ELSE BEGIN

    line          = ""
    READF, data_unit, line
    cols          = N_ELEMENTS(StrSplit(line))
    Point_Lun, data_unit, 0
    A             = DBLARR(cols, lines_per_layer)
    col_str       = STRTRIM(cols,2)
    fmt_str       = '(' + col_str + '(e24.16,1x),:)'

    SKIP_LUN, data_unit, l_start - 1, /LINES

    READF, data_unit, FORMAT = fmt_str, A

    A             = TRANSPOSE(A)

;   PRINT, "fetch_layer: size of A", SIZE(A)

;   PRINT, 'fetch_layer: n_lines         = ', n_lines
;   PRINT, 'fetch_layer: n3              = ', n3
;   PRINT, 'fetch_layer: cols            = ', cols    
;   PRINT, 'fetch_layer: lines_per_layer = ', lines_per_layer
;   PRINT, 'fetch_layer: l_start         = ', l_start
;   PRINT, 'fetch_layer: l_stop          = ', l_stop

    FREE_LUN, data_unit

;   FOR K = 0, lines_per_layer - 1 DO BEGIN

;    PRINT, FORMAT = '(a2,i4,a4,e24.16)', "A[", K, "] = ", A[K,1]
   
;   ENDFOR 

  ENDELSE


  RETURN, A

END