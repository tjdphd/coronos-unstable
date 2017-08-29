PRO ff_time_series, dsc_lab, fld_one

  i_time                     = 0  ; current time
  i_pe                       = 1  ; perpendicular kinetic energy
  i_ae                       = 2  ; perpendicular magnetic energy
  i_mo                       = 3  ; maximum vorticity
  i_imo                      = 4  ; location in z for maximum vorticity
  i_jmo                      = 5  ; location in plain for maximum vorticity
  i_oe                       = 6  ; square magnitude vorticity Laplacian
  i_mj                       = 7  ; maximum current
  i_imj                      = 8  ; location in z for maximum current
  i_jmj                      = 9  ; location in plain for maximum current
  i_ce                       = 10 ; square magnitude current Laplacian
  i_nuoe                     = 11 ; perpendicular viscous dissipation
  i_etace                    = 12 ; perpendicular resistive dissipation
  i_fe                       = 13 ; Poynting flux
  i_ftpow                    = 14 ; average injected foot-point energy per unit time
  i_ediss                    = 15 ; average cumulative dissipated energy per unit time
  i_deng                     = 16 ; average rage-of-change of total energy
  i_dedt                     = 17 ; instantaneous rate-of-change of total energy
  i_wtf                      = 18 ; average rate of energy gain minus energy loss
  i_consE                    = 19 ; consE
  i_tntauc                   = 20 ; current time in units of correlation time
  i_dt                       = 21 ; time increment
  i_dtvb                     = 22 ; dtvb
  i_ceocee                   = 23 ; Current density length scale 
  i_aveK2                    = 24 ; aveK2/t
  i_aveME                    = 25 ; time-averaged perpendicular magnetic field strength
  i_avepe                    = 26 ; time-averaged magnetic field footpoint velocity
  i_kze	                     = 27 ; parallel kinetic energy
  i_aze                      = 28 ; parallel magnetic energy
  i_ne                       = 29 ; internal energy
  i_bze                      = 30 ; Z Laplacian
  i_vze                      = 31 ; V_z Laplacian
  i_kapbze                   = 32 ; parallel conductive heat loss
  i_nuvze                    = 33 ; parallel viscous dissipation
   
  E_case                     = open_energy_data_file(dsc_lab)
  
  out_dev                    = 'X'
  
  ip1                        = scan_parameters('p1' , 0, desc_label ) ; power of 2 giving x-resolution
  ip2                        = scan_parameters('p2' , 0, desc_label ) ; power of 2 giving y-resolution
  n3                         = scan_parameters('p3' , 0, desc_label ) ; number of slices per data file
  mp                         = scan_parameters('np' , 0, desc_label ) ; number of processors used in run
  zl                         = scan_parameters('zl' , 0, desc_label ) ; total height along z of integration volume
   
  x_res                      = 2^ip1                                  ; resolution in x
  y_res                      = 2^ip2                                  ; resolution in y
  z_res                      = n3 * mp                                ; resolution in z
  
  i_x_res                    = UINT(x_res)
  i_y_res                    = UINT(y_res)
  i_z_res                    = UINT(z_res)
   
  str_x_res                  = STRTRIM(i_x_res, 2)
  str_y_res                  = STRTRIM(i_y_res, 2)
  str_z_res                  = STRTRIM(i_z_res, 2)
  
  str_zl                     = STRTRIM(UINT(zl),2)
  
  i_field                    = fld_one
  
   CASE i_field OF 
   '1': str_field            = 'pe'
   '2': str_field            = 'ae'
   '3': str_field            = 'mo'
   '4': str_field            = 'imo'
   '5': str_field            = 'jmo'
   '6': str_field            = 'oe'
   '7': str_field            = 'mj'
   '8': str_field            = 'imj'
   '9': str_field            = 'jmj'
  '10': str_field            = 'ce'
  '11': str_field            = 'nuoe'
  '12': str_field            = 'etace'
  '13': str_field            = 'fe'
  '14': str_field            = 'ftpow'
  '15': str_field            = 'ediss'
  '16': str_field            = 'deng'
  '17': str_field            = 'dedt'
  '18': str_field            = 'wtf'
  '19': str_field            = 'consE'
  '20': str_field            = 'tntauc'
  '21': str_field            = 'dt'
  '22': str_field            = 'dtvb'
  '23': str_field            = 'ceocee'
  '24': str_field            = 'aveK2'
  '25': str_field            = 'aveME'
  '26': str_field            = 'avepe'
  '27': str_field            = 'kze'
  '28': str_field            = 'aze'
  '29': str_field            = 'ne'
  '30': str_field            = 'bze'
  '31': str_field            = 'vze'
  '32': str_field            = 'kapbze'
  '33': str_field            = 'nuvze'
   ELSE: str_field           = 'Something_is_wrong'
   ENDCASE
  
  case_file_str              = 'ffts_' + str_field + '_' + str_x_res + '_' + str_z_res + '_zl_' + str_zl 

  res_str                    = str_x_res       + ' \times '   + str_y_res + ' \times '    + str_z_res
  case_res_str               = '$' + res_str   + ' (L = '     + str_zl    + ')'           + '$'

  E_size                     = SIZE(E_case)
  n_lines                    = E_size[1]
                             
  i_field                    = fld_one
                   
  max_y                      = MAX(E_case[*, i_field])
  min_y                      = MIN(E_case[*, i_field])
                   
  y_range                    = [min_y, max_y]
                             
  min_t                      = MIN(E_case[*, i_time])
  max_t                      = MAX(E_case[*, i_time])
                   
  t_range                    = [min_t, max_t]
                   
  cur_dir                    = GETENV('PWD')

  CASE i_field OF 
  '1': BEGIN
       str_title             = 'Total Perpendicular Kinetic Energy vs.t: '
       str_y_title           = '$E_{K,\perp}(t)$'
       str_plt_annote        = '$E_{K,\perp}(t) = (1/2)\int_0^L dz \int |\nabla_\perp\phi|^2 d\bf x\rm_\perp$'
       END
  '2': BEGIN
       str_title             = 'Total Perpendicular Magnetic Energy vs.t: '
       str_y_title           = '$E_{M,\perp}(t)$'
       str_plt_annote        = '$E_{M,\perp}(t) = (1/2)\int_0^L dz \int |\nabla_\perp A|^2 d\bf x\rm_\perp$'
       END
  '3': BEGIN
       str_title             = 'Maximum Vorticity vs.t: '
       str_y_title           = '$\Omega_{max}(t)$'
       str_plt_annote        = ''
       END
  '4': BEGIN
       str_title             = 'Should not be looking at this vs.t: '
       str_y_title           = '$forbidden(t)$'
       str_plt_annote        = ''
       END
  '5': BEGIN
       str_title             = 'Should not be looking at this vs.t: '
       str_y_title           = '$forbidden(t)$'
       str_plt_annote        = ''
       END
  '6': BEGIN
       str_title             = 'Integrated Square of Vorticity vs.t: '
       str_y_title           = '$I_\Omega$'
       str_plt_annote        = '$I_\Omega(t)=\int_0^L dz \int \Omega^2 d\bf x\rm_\perp$'
       END
  '7': BEGIN
       str_title             = 'Maximum Current vs.t: '
       str_y_title           = '$J_{max}(t)$'
       str_plt_annote        = ''
       END
  '8': BEGIN
       str_title             = 'Should not be looking at this vs.t: '
       str_y_title           = '$Forbidden(t)$'
       str_plt_annote        = ''
       END
  '9': BEGIN
       str_title             = 'Should not be looking at this vs.t: '
       str_y_title           = '$Forbidden(t)$'
       str_plt_annote        = ''
       END
 '10': BEGIN
       str_title             = 'Integrated Square of Current Density vs.t: '
       str_y_title           = '$I_J(t)$'
       str_plt_annote        = '$I_J(t) = \int_0^L dz \int J^2 d\bf x\rm_\perp$'
       END
 '11': BEGIN
       str_title             = 'Perpendicular Viscous Dissipation vs.t: '
       str_y_title           = '$W_{\nu,\perp}(t)$'
       str_plt_annote        = '$W_{\nu,\perp}(t)=\nu \int_0^L dz \int \Omega^2 d\bf x\rm_\perp$'
       END
 '12': BEGIN
       str_title             = 'Perpendicular Resistive Dissipation vs.t: '
       str_y_title           = '$W_{\eta,\perp}(t)$'
       str_plt_annote        = '$W_{\eta,\perp}(t) = \nu\int_0^L dz \int J^2 d\bf x\rm_\perp$'
       END
 '13': BEGIN
       str_title             = 'Exterior Boundary Poynting Flux vs.t: '
       str_y_title           = '$I(t)$'
       str_plt_annote        = '$I(t) =  \int\nabla_\perp A\cdot\nabla_\perp\phi d\bf x\rm_\perp|_{z=L}-' $
                                      + '\int\nabla_\perp A\cdot\nabla_\perp\phi d\bf x\rm_\perp|_{z=0}$'
       END
 '14': BEGIN
       str_title             = 'Average Injected Foot-Point Energy Per Unit Time vs.t: '
       str_y_title           = '$<I>(t)$'
       str_plt_annote        = '$<I>(t) = (\int_0^t I(t\prime) dt\prime)/t $'
       END
 '15': BEGIN
       str_title             = 'Time-Averaged Dissipated Energy Per Unit Time vs.t: '
       str_y_title           = '$<W>(t)$'
       str_plt_annote        = '$<W>(t) = [\int_0^t ( W_{\eta,\perp}(t\prime)+ W_{\nu,\perp}(t\prime) ) dt\prime]/t$'
       END
 '16': BEGIN
       str_title             = 'Average Rate-of-Change of Total Energy vs.t: '
       str_y_title           = '$<E>(t)$'
       str_plt_annote        = '$<E>(t)=[\int_0^t(dE_{K,\perp}/dt+dE_{M,\perp}/dt)dt\prime]/t$'
       END
 '17': BEGIN
       str_title             = 'Instantaneous Rate-of-Change of Total Energy vs.t: '
       str_y_title           = '$dE/dt(t)$'
       str_plt_annote        = '$dE/dt(t) = dE_{K,\perp}/dt +dE_{M,\perp}/dt$'
       END
 '18': BEGIN
       str_title             = 'Time-Average Energy Gains Minus Time-Average Energy Losses vs.t: '
       str_y_title           = '$<E(t)> - <W(t)>$'
       str_plt_annote        = ''
       END
 '19': BEGIN
       str_title             = 'consE vs.t: '
       str_y_title           = '$(dE/dt-dW/dt)dt$'
       str_plt_annote        = '$(dE/dt-dW/dt)dt=[I(t)-W_{\eta,\perp}(t)-W_{\nu,\perp}(t)-dE/dt]dt $'
       END
 '20': BEGIN
       str_title             = 'Current Time in Units of Correlation Time vs.t: '
       str_y_title           = '$t/\tau_C$'
       str_plt_annote        = ''
       END
 '21': BEGIN
       str_title             = 'dt vs.t: '
       str_y_title           = '$dt(t)$'
       str_plt_annote        = ''
       min_y                 = 0.0
       END
 '22': BEGIN
       str_title             = 'dtvb vs.t: '
       str_y_title           = '$dtvb(t)$'
       str_plt_annote        = ''
       END
 '23': BEGIN
       str_title             = 'Resistive Dissipation Length Scale vs.t: '
       str_y_title           = '$L_\eta(t)$'
       str_plt_annote        = '$L_\eta(t)=((\int_0^L dz \int J^2d\bf x\rm_\perp)/' $
                             + '(\int_0^Ldz\int|\nabla_\perp J|^2d\bf x\rm_\perp))^{1/2}$'
       END
 '24': BEGIN
       str_title             = 'Time-Averaged Inverse-Square, Current Density Length Scale vs.t: '
       str_y_title           = '$L^{-2}_J$'
       str_plt_annote        = '$L^{-2}_J = (\int_0^L dz \int J^2 d\bf x\rm_\perp)/' $
                             + '(\int_0^Ldz\int |\nabla_\perp A|^2d\bf x\rm_\perp)$'
       END
 '25': BEGIN
       str_title             = 'Time-Averaged Perpendicular Field Strength vs.t: '
       str_y_title           = '$<B_\perp>(t)$'
       str_plt_annote        = '$<B>(t)=(<E_{M,\perp}>(t))^{1/2}$' 
       END
 '26': BEGIN
       str_title             = 'Time-Averaged Photospheric Flow Velocity vs.t: '
       str_y_title           = '$<v_p>(t)$'
       str_plt_annote        = '$<v_p>(t)=(2(\int_0^t(\int |\nabla_\perp \phi|^2d\bf x\rm\perp)_{z=0} dt\prime)/t)^{1/2}$'
       END
 '27': BEGIN
       str_title             = 'Total Parallel Kinetic Energy vs.t: '
       str_y_title           = '$E_{K,\parallel}(t)$'
       str_plt_annote        = '$E_{K,\parallel}(t) = ?$'
       END
 '28': BEGIN
       str_title             = 'Total Parallel Magnetic Energy vs.t: '
       str_y_title           = '$E_{M,\parallel}(t)$'
       str_plt_annote        = '$E_{K,\parallel}(t) = ?$'
       END
 '29': BEGIN
       str_title             = 'Total Internal Energy vs.t: '
       str_y_title           = '$E_{th}(t)$'
       str_plt_annote        = '$E_{th}(t)=?$'
       END
 '30': BEGIN
       str_title             = 'Square Magnitude of Z vs.t: '
       str_y_title           = '$Z^2(t)$'
       str_plt_annote        = '$Z^2(t)=?$'
       END
 '31': BEGIN
       str_title             = 'Square Magnitude of V_z vs.t: '
       str_y_title           = '$V_z^2(t)$'
       str_plt_annote        = '$V_z^2(t)=?$'
       END
 '32': BEGIN
       str_title             = 'Parallel Conductive Heat Loss vs.t: '
       str_y_title           = '$E_{\kappa}(t)$'
       str_plt_annote        = '$E_{\kappa}(t)=?$'
       END
 '33': BEGIN
       str_title             = 'Parallel Viscous Dissipation vs.t: '
       str_y_title           = '$W_{\nu,\parallel}(t)$'
       str_plt_annote        = '$W_{\nu,\parallel}(t)=?$'
       END
  ELSE:BEGIN
        str_title            = 'Something_is_wrong'
       str_y_title           = '$E_X(t)$'
       str_plt_annote        = '$E_{K,\perp}(t) = (1\over 2)\int_0^Ldz \int |\nabla_\perp\phi|^2d{\bf x}_\perp$'
       END
  ENDCASE

  y_range                    = [min_y, max_y]

  max_ediss1                 = MAX(E_case[*,  i_ediss], i_max_ediss1)
  t_max_ediss1               =     E_case[i_max_ediss1, i_time]

  PRINT, '  max_ediss1 = ',   max_ediss1
  PRINT, 't_max_ediss1 = ', t_max_ediss1
  PRINT, 'i_max_ediss1 = ', i_max_ediss1

  AVL                        = E_case[*, i_field]

  i_beg                      = i_max_ediss1
  i_end                      = n_lines - 1

  IF (i_beg EQ i_end) THEN BEGIN
    i_beg                    = 0
  ENDIF

  ave                        = MEAN(  AVL[i_beg:i_end])
  sigma                      = STDDEV(AVL[i_beg:i_end])

  PRINT, 'average value for field:       ', i_field, ' is ', ave
  PRINT, 'standard deviation for field:  ', i_field, ' is ', sigma

  AVL[*]                     = ave

  str_title                  = str_title + ' For ' + case_res_str
  str_mean                   = STRTRIM(ave,2)
  len_mean                   = STRLEN(str_mean)

  eps_out = cur_dir + '/ffts/' + case_file_str + '.eps'

  X       = FLTARR(n_lines)
  Y       = FLTARR(n_lines)

  X       = E_case[0:(n_lines -1), i_time ]
  Y       = E_case[0:(n_lines -1), i_field]

; plt = PLOT( E_case[0:(n_lines-1),i_time], E_case[0:(n_lines-1),i_field],    $

; plt = PLOT( X, Y, /BUFFER)
; plt = PLOT( X, Y, /BUFFER,                                                  $

  plt    = PLOT( X, Y,                                                        $
                 AXIS_STYLE         = 2,                                      $
                 TITLE              = str_title,                              $
                 XTITLE             = "$t$",                                  $
                 YTITLE             = str_y_title,                            $
                 XRANGE             = t_range,                                $
                 YRANGE             = y_range,                                $
                 FONT_SIZE          = 14,                                     $
                 THICK              = 2                                       $
               )

 plt_ave = PLOT( E_case[0:(n_lines-1),i_time], AVL[0:(n_lines-1)],            $
                 /OVERPLOT,                                                   $
                 LINESTYLE          =   1,                                    $      
                 THICK              =   2,                                    $      
                 NAME               = "Mean = " + str_mean                    $      
               )
 
 
 plt_leg = LEGEND( TARGET           = [ plt_ave],                             $
;                  POSITION         = [0.5*t_range[1], 0.125*y_range[1]-0.01],     $
;                  POSITION         = [0.6*t_range[1], 0.725*y_range[1]],     $
                   POSITION         = [0.5*t_range[1], 0.5*y_range[1]],     $
                   /DATA,                                                     $
                   /AUTO_TEXT_COLOR                                           $
                 )

;plt_ann = TEXT( 0.025*t_range[1], 0.80*y_range[1], str_plt_annote,          $
 plt_ann = TEXT( 0.025*t_range[1], 0.80*y_range[1], str_plt_annote,           $
;plt_ann = TEXT( 0.0625*t_range[1], 0.20*y_range[1]-0.015, str_plt_annote,          $
;plt_ann = TEXT( 0.0625*t_range[1], 0.20*y_range[1], str_plt_annote,          $
                 /DATA,                                                       $
                 TARGET = [plt_ave],                                          $
                 FONT_SIZE = 16                                               $
               )

 plt.save, eps_out

END
