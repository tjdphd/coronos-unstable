PRO evsz, dsc_lab, efld, first_step, last_step, minmax_mode

  COMMON n_step, time

  !EXCEPT                 = 2

  ZERO                    = 0.0D                                      ; Define some constants for convenience
  ONE                     = 1.0D                                      ; everywhere
  TWO                     = 2.0D
  THREE                   = 3.0D
  PI                      = 3.141592653589793D
  TWO_THIRDS              = TWO / THREE
  TWO_PI                  = TWO*PI

  CZERO                   = DCOMPLEX(ZERO,ZERO)                       ; and to make sure that DOUBLE's are used
  CI                      = DCOMPLEX(ZERO,ONE)                        ; and to make sure that DOUBLE's are used

  i_all                   = 99
  i_z                     = 0 ; z - coordinate
  i_me                    = 1 ; megnetic energy in slab at z
  i_pe                    = 2 ; kinetic energy in slab at z
  i_ch                    = 3 ; cross helicity in slab at z
  i_ep                    = 4 ; elsasser energy
  i_em                    = 5 ; elsasser energy
  i_ce                    = 6 ; int (J^2 dz)
  i_oe                    = 7 ; int (Omega^2 dz)
  i_zp                    = 8 ; int( ( J + Omega)^2 dz)
  i_zm                    = 9 ; int( ( J - Omega)^2 dz)
  i_nch                   = 10 ; normalized cross helicity
  i_lm                    = 11 ; square - "magnetic" wave number 
  i_lp                    = 12 ; square - "kinetic"  wave number 
  i_lzp                   = 13 ; square - "Z+" wave number
  i_lzm                   = 14 ; square - "Z+" wave number
  i_te                    = 15 ; total energy at z
  i_re                    = 16 ; residual energy at z
  i_ne                    = 17 ; normalized residual energy
  i_rzp                   = 18 ; Z+ = |z+^2|^1/2  as determined from codes ep 
  i_rzm                   = 19 ; Z- = |z-^2|^1/2  as determined from codes em 
  i_zpm                   = 20 ; Z+Z- as determined from codes ep and em
  i_lpl                   = 21 ; length scale from F.T. of z^+
  i_lmi                   = 22 ; length scale from F.T. of z^+
  i_lpla                  = 23 ; alternative definition for lpl
  i_lmia                  = 24 ; alternative definition for lmi

; i_likp                  = 28 ; IK length length scale  lambda_IK,+  (see notes)
; i_likm                  = 29 ; IK length length scale  lambda_IK,-  (see notes)
; i_lkop                  = 30 ; Kolmogorov length scale lambda_Kol,+ (see notes)
; i_lkom                  = 31 ; Kolmogorov length scale lambda_Kol,- (see notes)

  i_likp                  = 25 ; IK length length scale  lambda_IK,+  (see notes)
  i_likm                  = 26 ; IK length length scale  lambda_IK,-  (see notes)
  i_lkop                  = 27 ; Kolmogorov length scale lambda_Kol,+ (see notes)
  i_lkom                  = 28 ; Kolmogorov length scale lambda_Kol,- (see notes)

; multi-threading required after here

  i_rzpf                  = 29 ; |z+^2|^1/2  as determined from ep calculated in post-processing
  i_rzmf                  = 30 ; |z-^2|^1/2  as determined from em calculated in post-processing
  i_zpmf                  = 31 ; Z+Z- as determined from ep and em calculated in post-processing

; i_rzpf                  = 25 ; |z+^2|^1/2  as determined from ep calculated in post-processing
; i_rzmf                  = 26 ; |z-^2|^1/2  as determined from em calculated in post-processing
; i_zpmf                  = 27 ; Z+Z- as determined from ep and em calculated in post-processing

  i_lsfp                  = 32 ; length scale from structure function for lambda +
  i_lsfm                  = 33 ; length scale from structure function for lambda -

  aumlaut                 = STRING(228B)
  CASE efld OF
  'zz'  : BEGIN 
          i_efld          = i_z
          str_efld        = "evsz_zz_"
          END
  'me'  : BEGIN 
          i_efld          = i_me
          str_efld        = "evsz_me_"
          str_title_pfx   = "Magnetic Energy vs z"
          str_y_title     = "$E_{M,\perp}(z)$"
          str_plt_annote  = "$E_{M,\perp}(z)=\int |\nabla_\perp A|^2 d\bf x\rm _\perp$"
          END
  'pe'  : BEGIN 
          i_efld          = i_pe
          str_efld        = "evsz_pe_"
          str_title_pfx   = "Kinetic Energy vs z"
          str_y_title     = "$E_{K,\perp}(z)$"
          str_plt_annote  = "$E_{K,\perp}(z)=\int |\nabla_\perp \phi|^2 d\bf x\rm _\perp$"
          END
  'ch'  : BEGIN 
          i_efld          = i_ch
          str_efld        = "evsz_ch_"
          str_title_pfx   = "Cross Helicity vs z"
          str_y_title     = "$H_C(z)$"
          str_plt_annote  = "$H_C(z)=2\int\nabla_\perp\phi\cdot\nabla_\perp A d\bf x\rm_\perp$"
          END
  'ep'  : BEGIN 
          i_efld          = i_ep
          str_efld        = "evsz_ep_"
          str_title_pfx   = "$Z^+$ Els" + aumlaut + "sser Energy vs z"
          str_y_title     = "$E_{Z}^+(z)$"
          str_plt_annote  = "$E_{Z}^+(z) = (1/2)\int|\nabla_\perp(\phi+A)|^2 d\bf x\rm_\perp $"
          END
  'em'  : BEGIN 
          i_efld          = i_em
          str_efld        = "evsz_em_"
          str_title_pfx   = "$Z^-$ Els" + aumlaut + "sser Energy vs z"
          str_y_title     = "$E_{Z}^{-}(z)$"
          str_plt_annote  = "$E_{Z}^{-}(z) = (1/2)\int|\nabla_\perp(\phi-A)|^2d\bf x\rm _\perp $"
          END
  'ce'  : BEGIN 
          i_efld          = i_ce
          str_efld        = "evsz_ce_"
          str_title_pfx   = "Square-Integrated Current vs z"
          str_y_title     = "$I_J(z)$"
          str_plt_annote  = "$I_J(z) = \int J^2 d\bf x\rm _\perp=\int |\nabla_\perp^2 A|^2d\bf x\rm_\perp$"
          END
  'oe'  : BEGIN 
          i_efld          = i_oe
          str_efld        = "evsz_oe_"
          str_title_pfx   = "Square-Integrated Vorticity vs z"
          str_y_title     = "$I_\Omega(z)$"
          str_plt_annote  = "$I_\Omega(z)=\int \Omega^2 d\bfx\rm_\perp=\int |\nabla_\perp^2\phi|^2d\bf x\rm_\perp$"
          END
  'zp'  : BEGIN 
          i_efld          = i_zp
          str_efld        = "evsz_zp_"
          str_title_pfx   = "Square-Integrated Els" +aumlaut + "sser Vorticity \omega_+ vs z"
          str_y_title     = "$I_{\omega^+}(z)$"
          str_plt_annote  = "$I_{\omega^+}(z)=\int \omega_+^2 d\bf x\rm_\perp =\int(\Omega+J)^2 d\bf x\rm_\perp $"
          END
  'zm'  : BEGIN 
          i_efld          = i_zm
          str_efld        = "evsz_zm_"
          str_title_pfx   = "Square-Integrated Els" +aumlaut+ " sser Vorticity \omega_- vs z"
          str_y_title     = "$I_{\omega^-}(z)$"
          str_plt_annote  = "$I_{\omega^-}(z)=\int \omega_-^2 d\bf x\rm_\perp = \int(\Omega -J)^2d\bf x\rm_\perp$"
          END
  'nch' : BEGIN 
          i_efld          = i_nch
          str_efld        = "evsz_nch_"
          str_title_pfx   = "Normalized Cross-Helicity vs z"
          str_y_title     = "$\sigma(z)$"
          str_plt_annote  = "$\sigma(z)=H_C(z)/(E_{K,\perp}(z)+E_{M,\perp}(z))$"
          END
  'lm'  : BEGIN 
          i_efld          = i_lm
          str_efld        = "evsz_lm_"
          str_title_pfx   = "Magnetic Energy Length Scale vs z"
          str_y_title     = "$L_M(z)$"
          str_plt_annote  = "$L_M(z)=|E_{M,\perp}(z)/I_{J}(z)|^{1/2}$"
          END
  'lp'  : BEGIN 
          i_efld          = i_lp
          str_efld        = "evsz_lp_"
          str_title_pfx   = "Kinetic Energy Length Scale vs z"
          str_y_title     = "$L_K(z)$"
          str_plt_annote  = "$L_K(z)=|E_{K,\perp}(z)/I_{\Omega}(z)|^{1/2}$"
          END
  'lzp' : BEGIN 
          i_efld          = i_lzp
          str_efld        = "evsz_lzp_"
          str_title_pfx   = "$Z^+$ Els" + aumlaut + "sser Length Scale vs z"
          str_y_title     = "$L_{Z^+}(z)$"
          str_plt_annote  = "$L_{Z^+}(z)=|E_{Z^+}(z)/I_{\omega^+}(z)|^{1/2}$"
          END
  'lzm' : BEGIN 
          i_efld          = i_lzm
          str_efld        = "evsz_lzm_"
          str_title_pfx   = "$Z^+$ Els" + aumlaut + "sser Length Scale vs z"
          str_y_title     = "$L_{Z^-}(z)$"
          str_plt_annote  = "$L_{Z^-}(z)=|E_{Z^-}(z)/J_{\omega^-}(z)|^{1/2}$"
          END
  'te'  : BEGIN 
          i_efld          = i_te
          str_efld        = "evsz_te_"
          str_title_pfx   = "Total Perpendicular Energy vs z"
          str_y_title     = "$E_{T,\perp}(z)$"
          str_plt_annote  = "$E_{T,\perp}(z)=E_{K,\perp}(z)+E_{M,\perp}(z)$"
          END
  're'  : BEGIN 
          i_efld          = i_re
          str_efld        = "evsz_re_"
          str_title_pfx   = "Residual Perpendicular Energy vs z"
          str_y_title     = "$E_{R,\perp}(z)$"
          str_plt_annote  = "$E_{R,\perp}(z)=E_{K,\perp}(z) - E_{M,\perp}(z)$"
          END
  'ne'  : BEGIN 
          i_efld          = i_ne
          str_efld        = "evsz_ne_"
          str_title_pfx   = "Normalized Residual Energy vs z"
          str_y_title     = "!7r!X!DD!N(z)"
          str_plt_annote  = "!7r!X!DD!N(z)=$E_{R,\perp}(z)/E_{T,\perp}(z)$"
          END
  'rzp' : BEGIN 
          i_efld          = i_rzp
          str_efld        = "evsz_rzp_"
          str_title_pfx   =  "$Z_+$ vs $z$"
          str_y_title     =  "$Z_+(z)$"
          str_plt_annote  =  "$Z_+=[E_{Z^+}(z)]^{1/2}$"
          END
  'rzm' : BEGIN 
          i_efld          = i_rzm
          str_efld        = "evsz_rzm_"
          str_title_pfx   = "$Z_-$ vs $z$"
          str_y_title     = "$Z_-(z)$"
          str_plt_annote  = "$Z_-(z)=[E_{Z^-}(z)]^{1/2}$"
          END
  'zpm' : BEGIN 
          i_efld          = i_zpm
          str_efld        = "evsz_zpm_"
          str_title_pfx   = "Z!D+!NZ!D-!N vs z"
          str_y_title     = "Z!D+!NZ!D-!N(z)"
          str_plt_annote  = "Z!D+!NZ!D-!N(z)=$[E_{Z^-}(z)E_{Z^+}]^{1/2}$"
          END
  'rzpf': BEGIN
          i_efld          = i_rzpf
          str_efld        = "evsz_rzp_pp"
          str_title_pfx   =  "$Z_+$ vs $z$ (determined in post-processing)"
          str_y_title     =  "$Z_+(z)$"
          str_plt_annote  =  "$Z_+=[E_{Z^+}(z)]^{1/2}$"
          END
  'rzmf': BEGIN
          i_efld          = i_rzmf
          str_efld        = "evsz_rzm_pp"
          str_title_pfx   = "$Z_-$ vs $z$ (determined in post-processing)"
          str_y_title     = "$Z_-(z)$"
          str_plt_annote  = "$Z_-(z)=[E_{Z^-}(z)]^{1/2}$"
          END
  'rzmf': BEGIN
          i_efld          = i_zpmf
          str_efld        = "evsz_zpm_pp"
          str_title_pfx   = "$Z_+Z_-$ vs $z$ (determined in post-processing) "
          str_y_title     = "$Z_+(z)Z_-(z)$"
          str_plt_annote  = "$Z_+(z)Z_-(z)=[E_{Z^-}(z)E_{Z^+}]^{1/2}$"
          END
  'likp': BEGIN 
          i_efld          = i_likp
          str_efld        = "evsz_likp_"
          str_title_pfx   = "IK Length Scale $\lambda^+_{IK}(z)$ vs $z$"
          str_y_title     = "$\lambda_{IK}(z)$"
          str_plt_annote  = "$\lambda_{IK}(z)=\pm(Z_+Z_-)^2/(dE_{Z^+}/dz)$"
          END
  'likm': BEGIN 
          i_efld          = i_likm
          str_efld        = "evsz_likm_"
          str_title_pfx   = "IK Length Scale $\lambda^-_{IK}(z)$ vs $z$"
          str_y_title     = "$L^-_{IK}(z)$"
          str_plt_annote  = "$L^-_{IK}(z)=\pm(Z_+Z_-)^2/(dE_{Z^-}/dz) $"
          END
  'lkop': BEGIN 
          i_efld          = i_lkop
          str_efld        = "evsz_lkop_"
          str_title_pfx   = "Kolmogorov Length Scale $\lambda^+_{Kol}(z)$ vs $z$"
          str_y_title     = "$\lambda^+_{Kol}(z)$"
          str_plt_annote  = "$\lambda^+_{Kol}(z)=\pm E_{Z^+}Z_-/(dE_{Z^+}/dz)$"
          END
  'lkom': BEGIN 
          i_efld          = i_lkom
          str_efld        = "evsz_lkom_"
          str_title_pfx   = "Kolmogorov Length Scale $\lambda^-_{Kol}(z)$ vs $z$"
          str_y_title     = "$\lambda^-_{Kol}(z)$"
          str_plt_annote  = "$\lambda^-_{Kol}(z)=\pm E_{Z^-}Z_+/(dE_{Z^-}/dz)$"
          END
  'lpl' : BEGIN 
          i_efld          = i_lpl
          str_efld        = "evsz_lpl_"
          str_title_pfx   = "Els"+aumlaut +"serr Length Scale $\lambda_+$ vs $z$"
          str_y_title     = "$\lambda_{+}(z)$"
          str_plt_annote  = "$\lambda_{+}(z)=(\pi/E_{Z^+}(z))[\Sigma_{\bf k\rm} |E_{Z^+}(\bf k\rm)|/|\bf k\rm|] $"
          END
  'lmi' : BEGIN 
          i_efld          = i_lmi
          str_efld        = "evsz_lmi_"
          str_title_pfx   =  "Els"+aumlaut +"serr Length Scale $\lambda_-$ vs $z$"
          str_y_title     = "$\lambda_{-}(z)$"
          str_plt_annote  = "$\lambda_{-}(z)=(\pi/E_{Z^-}(z))[\Sigma_{\bf k\rm} |E_{Z^-}(\bf k\rm)|/|\bf k\rm|] $"
          END
  'lpla': BEGIN 
          i_efld          = i_lpla
          str_efld        = "evsz_lpla_"
          str_title_pfx   = "Alternative Els" +aumlaut + "serr Length Scale $\lambda_{+,alt}$ vs $z$"
          str_y_title     = "$\lambda_{+,alt}$"
          str_plt_annote  = "$\lambda_{+,alt}=[(\pi/E_{Z^+}(z))\Sigma_{\bf k\rm} (|\phi(\bf k\rm)+A(\bf k\rm)|^2)]^{1/2}$"
          END
  'lmia': BEGIN 
          i_efld          = i_lmia
          str_efld        = "evsz_lmia_"
          str_title_pfx   =  "Alternative Els" +aumlaut + "serr Length Scale $\lambda_{-,alt}$ vs $z$"
          str_y_title     = "$\lambda_{-,alt}$"
          str_plt_annote  = "$\lambda_{-,alt}=[(\pi/E_{Z^-}(z))\Sigma_{\bf k\rm} (|\phi(\bf k\rm)-A(\bf k\rm)|^2)]^{1/2}$"
          END
  'lsfp': BEGIN 
          i_efld          = i_lsfp
          str_efld        = "evsz_lsfp_"
          str_title_pfx   = "Length Scale From Structure Function $\lambda_{s.f. +}$ vs $z$"
          str_y_title     = "$\lambda_{s.f. +}$"
          str_plt_annote  = "$\lambda_{s.f. +}=(\pi/Z^2_+(z))[\Sigma_{\bf k\rm} G(k^2L_\perp/2)|E_{Z^+}(\bf k\rm)|/|\bf k\rm|]$"
          END
  'lsfm': BEGIN 
          i_efld          = i_lsfm
          str_efld        = "evsz_lsfm_"
          str_title_pfx   = "Length Scale From Structure function $\lambda_{s.f. -}$ vs vs $z$"
          str_y_title     = "$\lambda_{s.f. -}(z)$"
          str_plt_annote  = "$\lambda_{s.f. -}(z)=(\pi/Z^2_-(z))[\Sigma_{\bf k\rm} (G(k^2L_\perp/2)|E_{Z^-}(\bf k\rm)|)/|\bf k\rm|]$"
          END
  'all' : BEGIN 
          i_efld          = i_all
          END
  ELSE: BEGIN 
        i_efld            = -1
        str_efld          = 'evsz_xx_' 
        str_title_pfx     = 'Something is terribly wrong'
        str_y_title       = 'This is just not right'
        str_plt_annote    = ''
        END
  ENDCASE

  cur_dir                 = GETENV('PWD')
  eps_out_dir             = cur_dir   + '/evsz/' + efld +    '/eps'
  
  ip1                     = scan_parameters('p1', 0, dsc_lab )                     ; power of 2 giving x-resolution
  ip2                     = scan_parameters('p2', 0, dsc_lab )                     ; power of 2 giving y-resolution
  n1                      = 2^ip1
  n2                      = 2^ip2
  n3                      = scan_parameters('p3' , 0, dsc_lab )                     ; number of slices per data file
  mp                      = scan_parameters('np' , 0, dsc_lab )                     ; number of processors used in run
  zl                      = scan_parameters('zl' , 0, dsc_lab )                     ; total height along z of integration volume
   
  x_res                   = n1                                                      ; resolution in x
  y_res                   = n2                                                      ; resolution in y
  z_res                   = n3 * mp                                                 ; resolution in z

  dz                      = zl / z_res
  dz_m1                   = z_res / zl                                              ; dz^(-1)

  cur_dir                 = GETENV('PWD')

  res_str                 = getResString(desc_label)
  gofx_dat                = cur_dir + '/evsz/' + 'gofx_' + res_str + '.dat'

  IF ( ( i_efld GE i_rzpf && i_efld LE i_lsfm) OR (i_efld EQ 99) ) THEN BEGIN ; need KK to be correct and maybe GOFX too

    IF ((i_efld GE i_lsfp) OR (i_efld EQ 99) ) THEN BEGIN                   ; need GOFX

;     res_str             = getResString(desc_label)
;     gofx_dat            = cur_dir + '/evsz/' + 'gofx_' + res_str + '.dat'

      IF (FILE_TEST(gofx_dat)) THEN BEGIN              ; If gofx  and kk are tabulated read them in
        
        PRINT, 'run_ae_esz: reading tabulated values for GOFX and KK'

        GOFX_LIN          = FLTARR(n1*n2)
        KK_LIN            = FLTARR(n1*n2)
        line              = FLTARR(2)

        OpenR, gofx_unit, gofx_dat, /GET_LUN
        FOR I = 0, n1*n2 - 1 DO BEGIN
          READF, gofx_unit, FORMAT = '(2(E24.16,1x),:/)', line
          GOFX_LIN[I]     = line[0]
          KK_LIN[I]       = line[1]
        ENDFOR

        GOFX              = REFORM(GOFX_LIN,n1,n2)
        KK                = REFORM(KK_LIN,n1,n2)
        
        GOFX_LIN          = REFORM(GOFX, n1*n2)

      ENDIF ELSE BEGIN                                   ; If gofx and kk are not tabulated then calculated them

        PRINT, 'run_ae_esz: calculating KK and GOFX'
        KK                = calcKK(n1,n2)
        GOFX              = calcGofx(KK, dsc_lab)

      ENDELSE
    ENDIF ELSE BEGIN                                     ; need KK but not GOFX

      KK                  = calcKK(n1,n2)                    ; calculate KK
      res_str             = getResString(desc_label)
;     gofx_dat            = cur_dir + '/ra_evsz/' + 'gofx_' + res_str + '.dat'
      gofx_dat            = cur_dir + '/evsz/' + 'gofx_' + res_str + '.dat'

      IF (FILE_TEST(gofx_dat)) THEN BEGIN                ; Go ahead read in GOFX and KK if they're available
        
        PRINT, 'run_ae_esz: reading tabulated values for GOFX and KK'

        GOFX_LIN          = FLTARR(n1*n2)
        KK_LIN            = FLTARR(n1*n2)
        line              = FLTARR(2)

        OpenR, gofx_unit, gofx_dat, /GET_LUN
        FOR I = 0, n1*n2 - 1 DO BEGIN
          READF, gofx_unit, FORMAT = '(2(E24.16,1x),:/)', line
          GOFX_LIN[I]     = line[0]
          KK_LIN[I]       = line[1]
        ENDFOR

        GOFX              = REFORM(GOFX_LIN,n1,n2)
        KK                = REFORM(KK_LIN,n1,n2)
        
        GOFX_LIN          = REFORM(GOFX, n1*n2)

      ENDIF ELSE BEGIN                                   ; Don't bother calculating GOFX if it's not tabulated
                                                         ; just give it a harmless definition
        GOFX              = FLTARR(n1,n2)
        GOFX[*,*]         = ONE

      ENDELSE
    ENDELSE
  ENDIF ELSE BEGIN
    IF (FILE_TEST(gofx_dat)) THEN BEGIN              ; If gofx  and kk are tabulated read them in
      
      PRINT, 'run_ae_esz: reading tabulated values for GOFX and KK'

      GOFX_LIN          = FLTARR(n1*n2)
      KK_LIN            = FLTARR(n1*n2)
      line              = FLTARR(2)

      OpenR, gofx_unit, gofx_dat, /GET_LUN
      FOR I = 0, n1*n2 - 1 DO BEGIN
        READF, gofx_unit, FORMAT = '(2(E24.16,1x),:/)', line
        GOFX_LIN[I]     = line[0]
        KK_LIN[I]       = line[1]
      ENDFOR

      GOFX              = REFORM(GOFX_LIN,n1,n2)
      KK                = REFORM(KK_LIN,n1,n2)
      
      GOFX_LIN          = REFORM(GOFX, n1*n2)

    ENDIF ELSE BEGIN                                   ; If gofx and kk are not tabulated then calculated them

      KK                = calcKK(n1,n2)                    ; calculate KK
      GOFX              = FLTARR(n1,n2)
      GOFX[*,*]         = ONE

    ENDELSE
  ENDELSE

  type_kk                 = TYPENAME(KK)

  IF (minmax_mode EQ  'global') THEN BEGIN
    
    glb_minmax_y          = FLTARR(i_lsfp,2)

    cur_dir               = GETENV('PWD')
    out_file              = cur_dir + '/evsz/' + 'evsz_global_extrema.out'
    OPENR, out_unit, out_file, /GET_LUN, ERROR = op_err

    IF ( (op_err EQ 0) AND (i_efld NE 99) ) THEN BEGIN
      ; read in extrema from file
      glb_minmax_y        = FLTARR(i_lsfm+1,2)
      next_line           = FLTARR(2)
      FOR K = 0, i_lsfm DO BEGIN
        READF, out_unit, FORMAT = '(2(e24.16,1x),:)', next_line
        glb_minmax_y[K,*] = next_line
      ENDFOR
      CLOSE, out_unit

    ENDIF ELSE BEGIN
;     glb_minmax_y             = global_efld_minmax(i_efld, dsc_lab, first_step, last_step, KK, GOFX)
      IF (i_efld NE 99) THEN BEGIN
        glb_minmax_y[i_efld,*] = global_efld_minmax(i_efld, dsc_lab, first_step, last_step, KK, GOFX)
;     ENDIF
      ENDIF ELSE BEGIN

;   ENDELSE
    
;   IF (i_efld EQ 99) THEN BEGIN

        glb_minmax_y             = global_efld_minmax(i_efld, dsc_lab, first_step, last_step, KK, GOFX)
        OPENW, out_unit, out_file, /GET_LUN, ERROR = op_err

        FOR K = 0, i_lsfm DO BEGIN

          PRINTF, out_unit, FORMAT = '(2(e24.16,1x),:)', glb_minmax_y[K,*]

        ENDFOR
        CLOSE, out_unit

        READ, Q, PROMPT="completed determination of global extrema. Select 0 to quit or select new field to plot:"
        IF ( Q EQ 0) THEN BEGIN
          EXIT
        ENDIF ELSE BEGIN

          i_efld = Q
          PRINT, "evsz: you selected i_efld = ", i_efld

        ENDELSE

      ENDELSE
    ENDELSE

;   ENDIF

  ENDIF

  FOR I = first_step, last_step DO BEGIN

    E                     = open_evsz_data_file( dsc_lab, I)
    extE                  = calc_extE(E, dsc_lab, I, i_efld, KK, GOFX )

    ineq                  = WHERE(E[*,0:i_lmia] NE extE[*,0:i_lmia], count)

    IF (count NE 0) THEN BEGIN
      PRINT, ' evsz: WARNING - E and extE inconsistent. count = ', count, ' at time ', time
      STOP
    ENDIF

    E                     = extE
    size_E                = SIZE(extE)
    n_lines               = size_E[1]
    n_cols                = size_E[2]

;   dzpdz                 = FLTARR(n_lines)
;   dzmdz                 = FLTARR(n_lines)

;   FOR J = 0, n_lines - 2 DO BEGIN

;     delta_z             = E[J+1,i_z] - E[J, i_z]

;     IF (delta_z NE dz) THEN BEGIN
;       PRINT, 'ra_evsz: WARNING - value of dz inconsistent with deltaE[*,0] = ',  delta_z
;     ENDIF

;     dzpdz[J]            = (E[J+1,i_ep] - E[J, i_ep]) * dz_m1
;     dzmdz[J]            = (E[J+1,i_em] - E[J, i_em]) * dz_m1

;   ENDFOR

;   dzpdz[n_lines - 1]    = dzpdz[n_lines - 2]
;   dzmdz[n_lines - 1]    = dzmdz[n_lines - 2]

;   pzeros                = WHERE( dzpdz EQ ZERO, pzcount )

    ; Kolmogorov: 

;   E[*,i_lkop]           = ABS((E[*, i_ep] * E[*,i_rzm]) / dzpdz[*])
;   E[*,i_lkom]           = ABS((E[*, i_em] * E[*,i_rzp]) / dzmdz[*])


    ; (ep^1/2 * em^1/2)

;   E[*, i_zpm]           = E[*, i_rzp] * E[*, i_rzm]

    ; Iroshnikov-Kraichnan

;   E[*, i_likp]          = ABS(E[*,i_zpm]^2 / dzpdz[*])
;   E[*, i_likm]          = ABS(E[*,i_zpm]^2 / dzmdz[*])

    i_x_res               = UINT(x_res)
    i_y_res               = UINT(y_res)
    i_z_res               = UINT(z_res)
     
    str_x_res             = STRTRIM(i_x_res, 2)
    str_y_res             = STRTRIM(i_y_res, 2)
    str_z_res             = STRTRIM(i_z_res, 2)
    str_zl                = STRTRIM(UINT(zl),2)
    case_res_str          = ' (L = ' + str_zl    + ')'

    str_time              = STRTRIM(time,2)

    str_step              = STRTRIM(I,2)
    str_last_step         = STRTRIM(last_step,2)
    max_z_digits          = STRLEN(str_last_step)
    IF (max_z_digits LT 3) THEN max_z_digits = 3
    z_digits              = STRLEN(str_step)
    zero_pad              = max_z_digits - z_digits
    zero_str              = ''
    IF(zero_pad NE 0) THEN BEGIN
      FOR J = 1, zero_pad DO BEGIN
        zero_str          = zero_str + '0'
      ENDFOR
    ENDIF

    str_step              = zero_str + str_step

    case_file_str         = str_efld + str_x_res + '_' + str_z_res + '_zl_' + str_zl + '-' + str_step
    eps_out               = eps_out_dir + '/' + case_file_str + '.eps'
    tif_out               = eps_out_dir + '/' + case_file_str + '.tif'

    str_title             = str_title_pfx + ' at t = ' + str_time + case_res_str

    min_x                 = MIN(E[*, i_z   ])
    max_x                 = MAX(E[*, i_z   ])

    IF (min_x GT 0.0) THEN min_x = 0.0

    x_rng                 = [min_x, max_x]

    IF (minmax_mode EQ 'global') THEN BEGIN
      y_rng               = glb_minmax_y[i_efld,*]
    ENDIF ELSE BEGIN
      min_y               = MIN(E[*, i_efld])
      max_y               = MAX(E[*, i_efld])
      y_rng               = [min_y, max_y]
    ENDELSE

    X                     = FLTARR(n_lines)
    Y                     = FLTARR(n_lines)

    X[*]                  = E[*, i_z   ]
    Y[*]                  = E[*, i_efld]

    plt                   = PLOT( X, Y,                                       $
                                 /BUFFER,                                     $
                                 AXIS_STYLE = 2,                              $
                                 TITLE      = str_title,                      $
                                 XTITLE     = 'z',                            $
                                 YTITLE     = str_y_title,                    $
                                 XRANGE     = x_rng,                          $
                                 YRANGE     = y_rng,                          $
                                 FONT_SIZE  = 14                              $
                                )



    plt_ann               = TEXT(0.025*x_rng[1],0.75*y_rng[1],str_plt_annote, $
                                 /DATA,                                       $
                                 TARGET     = [plt],                          $
                                 FONT_SIZE  = 16                              $
                                )
    plt.save, eps_out
    plt.erase

  ENDFOR

END
