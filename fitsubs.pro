pro Ha, x, par, f, pder

z=par[0]
ww-par[1]

Ha = 6562.8*(1+z)
NII_b = 6548.1*(1+z)
NII_r = 6583*(1+z)
Hainten = par[2]
NIIinten = par[3]
cont = par[4]

Ha_ln = (Hainten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-Ha)^2/ww^2))
NII_a_ln = (NIIinten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-NII_r)^2/ww^2))
NII_b_ln = (NIIinten/3./ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-NII_b)^2/ww^2))

f = cont + Ha_ln + NII_a_ln + NII_b_ln

pder = fltarr(n_elements(x),n_elements(par))   ; no value returned.

end


;---------------------------------------------------------------------------------------
pro Ha_br, x, par, f, pder

z=par[0]
ww-par[1]

Ha = 6562.8*(1+z)
NII_b = 6548.1*(1+z)
NII_r = 6583*(1+z)
Hainten = par[2]
NIIinten = par[3]
cont = par[4]

z_br = par[5]
ww_br = par[6]
Ha_br = 6562.8*(1+z_br)
NII_b_br = 6548.1*(1+z_br)
NII_r_br = 6583*(1+z_br)
Hainten_br = par[7]
NIIinten_br = par[8]



Ha_ln = (Hainten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-Ha)^2/ww^2))
NII_a_ln = (NIIinten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-NII_r)^2/ww^2))
NII_b_ln = (NIIinten/3./ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-NII_b)^2/ww^2))

Ha_br_ln = (Hainten_br/ww_br/sqrt(2.0*!pi)) * (exp(-0.5*(x-Ha_br)^2/ww^2))
NII_a_br_ln = (NIIinten_br/ww_br/sqrt(2.0*!pi)) * (exp(-0.5*(x-NII_r_br)^2/ww^2))
NII_b_br_ln = (NIIinten_br/3./ww_br/sqrt(2.0*!pi)) * (exp(-0.5*(x-NII_b_br)^2/ww^2))

f = cont + Ha_ln + NII_a_ln + NII_b_ln + Ha_br_ln + NII_a_br_ln + NII_b_br_ln

pder = fltarr(n_elements(x),n_elements(par))   ; no value returned.

end


;---------------------------------------------------------------------------------------
pro Hb, x, par, f, pder

z=par[0]
ww-par[1]

Hb = 4861.32*(1+z)
OIII_b = 4959.00*(1+z)
OIII_r = 5007.00*(1+z)
Hbinten = par[2]
OIIIinten = par[3]
cont = par[4]

Hb_ln = (Hbinten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-Hb)^2/ww^2))
OIII_a_ln = (OIIIinten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-OIII_r)^2/ww^2))
OIII_b_ln = (OIIIinten/3./ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-OIII_b)^2/ww^2))

f = cont + Hb_ln + OIII_a_ln + OIII_b_ln

pder = fltarr(n_elements(x),n_elements(par))   ; no value returned.

end


;---------------------------------------------------------------------------------------
pro Hb_br, x, par, f, pder

z=par[0]
ww-par[1]

Hb = 4861.32*(1+z)
OIII_b = 4959.00*(1+z)
OIII_r = 5007.00*(1+z)
Hbinten = par[2]
OIIIinten = par[3]
cont = par[4]

z_br = par[5]
ww_br = par[6]
Hb_br = 4861.32*(1+z_br)
OIII_b_br = 4959.00*(1+z_br)
OIII_r_br = 5007.00*(1+z_br)
Hbinten_br = par[2]
OIIIinten_br = par[3]



Hb_ln = (Hbinten_br/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-Hb)^2/ww^2))
OIII_a_ln = (OIIIinten_br/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-OIII_r)^2/ww^2))
OIII_b_ln = (OIIIinten_br/3./ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-OIII_b)^2/ww^2))

Hb_ln_br = (Hbinten/ww_br/sqrt(2.0*!pi)) * (exp(-0.5*(x-Hb_br)^2/ww_br^2))
OIII_a_ln_br = (OIIIinten/ww_br/sqrt(2.0*!pi)) * (exp(-0.5*(x-OIII_r_br)^2/ww_br^2))
OIII_b_ln_br = (OIIIinten/3./ww_br/sqrt(2.0*!pi)) * (exp(-0.5*(x-OIII_b_br)^2/ww_br^2))

f = cont + Hb_ln + OIII_a_ln + OIII_b_ln + Hb_ln_br + OIII_a_ln_br + OIII_b_ln_br

pder = fltarr(n_elements(x),n_elements(par))   ; no value returned.

end


;---------------------------------------------------------------------------------------
pro OI, x, par, f, pder

z=par[0]
ww-par[1]

OI_b = 6300.00*(1+z)
OI_r = 6366.00*(1+z)
OI_b_inten = par[2]
OI_r_inten = par[3]
cont = par[4]

OI_b_ln = (OI_b_inten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-OI_b)^2/ww^2))
OI_r_ln = (OI_r_inten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-OI_r)^2/ww^2))

f = cont + OI_b_ln + OI_r_ln

pder = fltarr(n_elements(x),n_elements(par))   ; no value returned.

end


;---------------------------------------------------------------------------------------
pro OI_br, x, par, f, pder

z=par[0]
ww-par[1]

OI_b = 6300.00*(1+z)
OI_r = 6366.00*(1+z)
OI_b_inten = par[2]
OI_r_inten = par[3]
cont = par[4]

z_br = par[5]
ww_br = par[6]
OI_b_br = 6300.00*(1+z_br)
OI_r_br = 6366.00*(1+z_br)
OI_b_inten_br = par[2]
OI_r_inten_br = par[3]


OI_b_ln = (OI_b_inten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-OI_b)^2/ww^2))
OI_r_ln = (OI_r_inten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-OI_r)^2/ww^2))

OI_b_br_ln = (OI_b_inten_br/ww_br/sqrt(2.0*!pi)) * (exp(-0.5*(x-OI_b_br)^2/ww_br^2))
OI_r_br_ln = (OI_r_inten_br/ww_br/sqrt(2.0*!pi)) * (exp(-0.5*(x-OI_r_br)^2/ww_br^2))

f = cont + OI_b_ln + OI_r_ln + OI_b_ln_br + OI_r_ln_br

pder = fltarr(n_elements(x),n_elements(par))   ; no value returned.

end

;---------------------------------------------------------------------------------------
pro SII, x, par, f, pder

z=par[0]
ww-par[1]

SII_b = 6716.00*(1+z)
SII_r = 6731.00*(1+z)
SII_b_inten = par[2]
SII_r_inten = par[3]
cont = par[4]

SII_b_ln = (SII_b_inten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-SII_b)^2/ww^2))
SII_r_ln = (SII_r_inten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-SII_r)^2/ww^2))

f = cont + SII_b_ln + SII_r_ln

pder = fltarr(n_elements(x),n_elements(par))   ; no value returned.

end


;---------------------------------------------------------------------------------------
pro SII_br, x, par, f, pder

z=par[0]
ww-par[1]

SII_b = 6716.00*(1+z)
SII_r = 6731.00*(1+z)
SII_b_inten = par[2]
SII_r_inten = par[3]
cont = par[4]

z_br = par[5]
ww_br = par[6]
SII_b_br = 6716.00*(1+z_br)
SII_r_br = 6731.00*(1+z_br)
SII_b_inten_br = par[2]
SII_r_inten_br = par[3]


SII_b_ln = (SII_b_inten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-SII_b)^2/ww^2))
SII_r_ln = (SII_r_inten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-SII_r)^2/ww^2))

SII_b_br_ln = (SII_b_inten_br/ww_br/sqrt(2.0*!pi)) * (exp(-0.5*(x-SII_b_br)^2/ww_br^2))
SII_r_br_ln = (SII_r_inten_br/ww_br/sqrt(2.0*!pi)) * (exp(-0.5*(x-SII_r_br)^2/ww_br^2))

f = cont + SII_b_ln + SII_r_ln + SII_b_ln_br + SII_r_ln_br

pder = fltarr(n_elements(x),n_elements(par))   ; no value returned.

end




