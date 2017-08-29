FUNCTION getSfldIndex, efld


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

; E extensions starts here

  i_likp                  = 25 ; IK length length scale  lambda_IK,+  (see notes)
  i_likm                  = 26 ; IK length length scale  lambda_IK,-  (see notes)
  i_lkop                  = 27 ; Kolmogorov length scale lambda_Kol,+ (see notes)
  i_lkom                  = 28 ; Kolmogorov length scale lambda_Kol,- (see notes)

; multi-threading required after here

  i_rzpf                  = 29 ; |z+^2|^1/2  as determined from ep calculated in post-processing
  i_rzmf                  = 30 ; |z-^2|^1/2  as determined from em calculated in post-processing
  i_zpmf                  = 31 ; Z+Z- as determined from ep and em calculated in post-processing
  i_lsfp                  = 32 ; length scale from structure function for lambda +
  i_lsfm                  = 33 ; length scale from structure function for lambda -
  i_all                   = 99

  CASE efld OF
  'zz'  : BEGIN 
          i_efld          = i_z
  'me'  : BEGIN 
          i_efld          = i_me
          END
  'pe'  : BEGIN 
          i_efld          = i_pe
          END
  'ch'  : BEGIN 
          i_efld          = i_ch
          END
  'ep'  : BEGIN 
          i_efld          = i_ep
          END
  'em'  : BEGIN 
          i_efld          = i_em
          END
  'ce'  : BEGIN 
          i_efld          = i_ce
          END
  'oe'  : BEGIN 
          i_efld          = i_oe
          END
  'zp'  : BEGIN 
          i_efld          = i_zp
          END
  'zm'  : BEGIN 
          i_efld          = i_zm
          END
  'nch' : BEGIN 
          i_efld          = i_nch
          END
  'lm'  : BEGIN 
          i_efld          = i_lm
          END
  'lp'  : BEGIN 
          i_efld          = i_lp
          END
  'lzp' : BEGIN 
          i_efld          = i_lzp
          END
  'lzm' : BEGIN 
          i_efld          = i_lzm
          END
  'te'  : BEGIN 
          i_efld          = i_te
          END
  're'  : BEGIN 
          i_efld          = i_re
          END
  'ne'  : BEGIN 
          i_efld          = i_ne
          END
  'rzp' : BEGIN 
          i_efld          = i_rzp
          END
  'rzm' : BEGIN 
          i_efld          = i_rzm
          END
  'zpm' : BEGIN 
          i_efld          = i_zpm
          END
  'rzpf': BEGIN
          i_efld          = i_rzpf
          END
  'rzmf': BEGIN
          i_efld          = i_rzmf
          END
  'rzmf': BEGIN
          i_efld          = i_zpmf
          END
  'likp': BEGIN 
          i_efld          = i_likp
          END
  'likm': BEGIN 
          i_efld          = i_likm
          END
  'lkop': BEGIN 
          i_efld          = i_lkop
          END
  'lkom': BEGIN 
          i_efld          = i_lkom
          END
  'lpl' : BEGIN 
          i_efld          = i_lpl
          END
  'lmi' : BEGIN 
          i_efld          = i_lmi
          END
  'lpla': BEGIN 
          i_efld          = i_lpla
          END
  'lmia': BEGIN 
          i_efld          = i_lmia
          END
  'lsfp': BEGIN 
          i_efld          = i_lsfp
          END
  'lsfm': BEGIN 
          i_efld          = i_lsfm
          END
  'all' : BEGIN 
          i_efld          = i_all
          END
  ELSE  : BEGIN 
          i_efld            = -1
          END

  ENDCASE

  RETURN, i_efld
END
