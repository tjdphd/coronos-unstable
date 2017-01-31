FUNCTION open_energy_data_file, desc_label, nf

prefix       = scan_parameters('prefix',   0, desc_label )
ip1          = scan_parameters('p1',       0, desc_label )                     ; power of 2 giving x-resolution
ip2          = scan_parameters('p2',       0, desc_label )                     ; power of 2 giving y-resolution
n3           = scan_parameters('p3' ,      0, desc_label )                     ; number of slices per data file
mp           = scan_parameters('np' ,      0, desc_label )                     ; number of processors used in run
zl           = scan_parameters('zl' ,      0, desc_label )                     ; total height along z of integration volume
data_dir     = scan_parameters('data_dir', 0, desc_label )
 
x_res        = 2^ip1                                                      ; resolution in x
y_res        = 2^ip2                                                      ; resolution in y
z_res        = n3 * mp                                                    ; resolution in z

i_x_res      = UINT(x_res)
i_y_res      = UINT(y_res)
i_z_res      = UINT(z_res)
 
str_x_res    = STRTRIM(i_x_res, 2)
str_z_res    = STRTRIM(i_z_res, 2)

str_res_lab  = str_x_res + '_' + str_z_res

data_file    = prefix +'_' + str_res_lab + '.' + 'o' + desc_label            ; 'rmct2_64_16.ots'

cur_dir      = GETENV('PWD')
data_file    = cur_dir + '/' + data_dir + '/' + data_file

PRINT, 'opening file: ', data_file
OPENR, data_unit, data_file, /GET_LUN

nlines       = FILE_LINES(data_file)
line         = ""
READF, data_unit, line
cols         = N_ELEMENTS(StrSplit(line))
str_cols     = STRING(cols)
str_fmt      = '(' + str_cols + '(e24.20,1x),:)'
Point_Lun, data_unit, 0

Eline        = FLTARR(cols, 1)
E            = FLTARR(cols, nlines)

FOR I = 0, nlines - 1 DO BEGIN

   READF, data_unit, FORMAT = str_fmt, Eline
   E[*, I]   = Eline
  
ENDFOR

 E           = TRANSPOSE(E)  

FREE_LUN, data_unit
RETURN, E

END
