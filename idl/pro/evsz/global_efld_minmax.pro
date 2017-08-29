FUNCTION global_efld_minmax, i_efld, dsc_lab, first_step, last_step, KK, GOFX

  total_steps                               = last_step - first_step + 1

  IF (i_efld NE 99) THEN BEGIN
    efld_minmax_by_step                     = FLTARR(2, total_steps)

    FOR I = first_step, last_step DO BEGIN

      E                                     = open_evsz_data_file( dsc_lab, I)
      extE                                  = calc_extE( E, dsc_lab , I, i_efld, KK, GOFX)

      min_E_by_step                         = MIN(extE[*,i_efld])
      max_E_by_step                         = MAX(extE[*,i_efld])

      efld_minmax_by_step[0,I - first_step] = min_E_by_step 
      efld_minmax_by_step[1,I - first_step] = max_E_by_step 

    ENDFOR

    glb_efld_min                            = MIN(efld_minmax_by_step[0,*])
    glb_efld_max                            = MAX(efld_minmax_by_step[1,*])

    glb_minmax                              = [glb_efld_min, glb_efld_max]
 
  ENDIF ELSE BEGIN

    PRINT, "global_efld_minmax: calculating extrema for all fields..."

    efld_minmax_by_step                     = FLTARR(34, 2, total_steps )
    min_E_by_step                           = FLTARR(34,    total_steps ) 
    max_E_by_step                           = FLTARR(34,    total_steps )
    glb_minmax                              = FLTARR(34, 2              )

    FOR I = first_step, last_step DO BEGIN

      E                                     = open_evsz_data_file( dsc_lab, I)
      extE                                  = calc_extE( E, dsc_lab , I, i_efld, KK, GOFX)

      FOR J = 0, 33 DO BEGIN

        min_E_by_step[J,I-first_step]     = MIN( extE[*,J] )
        max_E_by_step[J,I-first_step]     = MAX( extE[*,J] )

      ENDFOR

      efld_minmax_by_step[*, 0, I-first_step] = min_E_by_step[*,I-first_step]
      efld_minmax_by_step[*, 1, I-first_step] = max_E_by_step[*,I-first_step]

    ENDFOR
    
    FOR J = 0, 33 DO BEGIN

      glb_minmax[J, 0]                    = MIN(efld_minmax_by_step[J, 0, *])
      glb_minmax[J, 1]                    = MAX(efld_minmax_by_step[J, 1, *])

    ENDFOR

  ENDELSE

  RETURN, glb_minmax

END
