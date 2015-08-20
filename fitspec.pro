; ---------------------------------------------------------------------------------------------------------------

pro fitspec

fits_read,'coadd-HaNII.fits',finalcube
restore,'wl_HaNII.dat'
fits_read, 'sky_HaNII.fits', sky


sig2d = finalcube[*,*,373:613]
finalcube = finalcube[*,*,100:340]
wl = wl[100:340]
sz = size(finalcube)
sky = sky[100:340]
var = 1/sky

fitcube = fltarr(sz[1],sz[2],sz[3])     ; holds best fit
acube = fltarr(sz[1],sz[2],5)           ; holds best parameters
bcube = fltarr(sz[1],sz[2],3)           ; holds best broad line parameters
aerrorcube = fltarr(sz[1],sz[2],8)      ; holds +- 1 sigma values for vel.
;askycube = fltarr(sz[1],sz[2])          ; holds position of the sky
;line.
alast = [0.034354,0.1,2.0,0,0.1] ; [z, Ha inten, Ha width, NII inten, contin]
brlf =  [0.034354,0.5,2.0,0,0.3, 1,10,0.034354] ; alast plus [bl inten, bl width, bl z]

; looping in this order makes the fit wander less between columns.

;window, 0, xsize=700,ysize=700, retain=2
;setcolours
!p.multi=[0,3,3]
;  just loop over a smll region for now.
akeep = alast
bkeep = brlf
for j=3,sz[2]-3 do for i = 3,sz[1]-3 do begin
   a = alast
   b = brlf
   if j eq 1 then a=akeep
   if j eq 1 then b=bkeep

i=35
j=45

;  j=20
;  i=20

;  i=10
;  j=10

  x = wl
  y = fltarr(sz[3])    			; make this average of 9 lenslets.
  sig = fltarr(sz[3])
  for k = 0,sz[3]-1 do y(k) = total(finalcube[i,j,k])
  for k = 0,sz[3]-1 do sig(k) = total(sig2d[i,j,k])

  sig2 = stddev(sig)
  sig2 = sig2*sig2

  w = fltarr(sz[3])+1.
  ;formpfit:  w = sky/median(sky)
  w = var/median(var) ; bodge to get chi^2=1 in noise


  a2 = a

  av =    curvefit(x, y, w, a2, sigma, chisq=bestnorme1,function_name="func",/noderivative)    ;(moment(y))[0]


  chi0 = total((w*(y - av)^2)/sig2)
  chi0r = chi0/n_elements(y)
  print, 'Total chi^2 when no emission line = ', chi0
  print, 'Reduced chi^2 when no emission line = ', chi0r

;stop

  fit = curvefit(x,y,w,a,sigma,chisq=bestnorme2,function_name="triplet",/noderivative)
  chi1 = total((w*(y - fit)^2)/sig2)
  chi1r = chi1/(n_elements(y)-(n_elements(a)-1))
  print, 'Total chi^2 when fit inclubed = ', chi1
  print, 'Reduced chi^2 when fit included = ', chi1r

;  stop

;  plot,x,y
;  oplot,x,fit,color=2

  chi2lim = 49.


;compute p
if ( alpha gt p ) then begin


  if (chi0-chi1 gt chi2lim) then begin
    acube[i,j,*] = a
  endif else   acube[i,j,*] = a2


endif else   acube[i,j,*] = a2

;to remocve
  if (chi0-chi1 gt chi2lim) then begin
     ;plot, wl, y
     ;oplot, wl, fit, color=2
     ;fitcube[i,j,*] = fit
     ;acube[i,j,*] = a
     ;print, 'best parameters = ', a
     ;print, '  with ', chi0-chi1, ' at ', i, j


     fitb = curvefit(x,y,w,b,sigma,chisq=chi2,function_name="tripletbr",/noderivative)
     chi2 = total((w*(y - fitb)^2)/sig2)
     chi2r = chi1/(n_elements(y)-(n_elements(b)-1))
     print, 'Total chi^2 when broad line inclubed = ', chi2
     print, 'Reduced chi^2 when broad line included = ', chi2r

     lim = chi1 - chi2

     if lim ge 64 then begin
       plot, wl, y
       oplot, wl, fitb, color=2
       fitcube[i,j,*] = fitb
       acube[i,j,*] = b[0:4]
       bcube[i,j,*] = b[5:7]
       print, 'best parameters = ', b
       print, '  with ', chi1-chi2, ' at ', i, j

      ;err1=errsinglet(wl,y,b,w,chi2,sig2,[0.00001,0,0,0,0,0,0,0],0,brtriplet=1)
      ;err2=errsinglet(wl,y,b,w,chi2,sig2,[0,0.001,0,0,0,0,0,0],1,brtriplet=1)
     ;err3=errsinglet(wl,y,b,w,chi2,sig2,[0,0,0.00001,0,0,0,0,0],2,brtriplet=1)
     ;err4=errsinglet(wl,y,b,w,chi2,sig2,[0,0,0,0.00001,0,0,0,0],3,brtriplet=1)
     ;err5=errsinglet(wl,y,b,w,chi2,sig2,[0,0,0,0,0.00001,0,0,0],4,brtriplet=1)
     ;err6=errsinglet(wl,y,b,w,chi2,sig2,[0,0,0,0,0,0.001,0,0],5,brtriplet=1)
     ;err7=errsinglet(wl,y,b,w,chi2,sig2,[0,0,0,0,0,0,0.0001,0],6,brtriplet=1)
     ;err8=errsinglet(wl,y,b,w,chi2,sig2,[0,0,0,0,0,0,0,00.00001],7,brtriplet=1)

     ;aerrorcube[i,j,0] = err1
     ;aerrorcube[i,j,1] = err2
     ;aerrorcube[i,j,2] = err3
     ;aerrorcube[i,j,3] = err4
     ;aerrorcube[i,j,4] = err5
     ;aerrorcube[i,j,5] = err6
     ;aerrorcube[i,j,6] = err7
     ;aerrorcube[i,j,7] = err8

     endif else begin

         plot, wl, y
         oplot, wl, fit, color=2
         oplot, wl, fitb, color=4
         fitcube[i,j,*] = fit
         acube[i,j,*] = a
         print, 'best parameters = ', a
         print, '  with ', chi0-chi1, ' at ', i, j

      ;err1=errsinglet(wl,y,a,w,chi1,sig2,[0.00001,0,0,0,0],0,triplet=1)
      ;err2=errsinglet(wl,y,a,w,chi1,sig2,[0,0.001,0,0,0],1,triplet=1)
     ;err3=errsinglet(wl,y,a,w,chi1,sig2,[0,0,0.00001,0,0],2,triplet=1)
     ;err4=errsinglet(wl,y,a,w,chi1,sig2,[0,0,0,0.00001,0],3,triplet=1)
     ;err5=errsinglet(wl,y,a,w,chi1,sig2,[0,0,0,0,0.00001],4,triplet=1)


     ;aerrorcube[i,j,0] = err1
     ;aerrorcube[i,j,1] = err2
     ;aerrorcube[i,j,2] = err3
     ;aerrorcube[i,j,3] = err4
     ;aerrorcube[i,j,4] = err5
     ;aerrorcube[i,j,5] = 0
     ;aerrorcube[i,j,6] = 0
     ;aerrorcube[i,j,7] = 0

     endelse



;     stop

;wait,5

  endif  else  begin
    for k = 0,sz[3]-1 do y(k) = total(finalcube[i-1:i+1,j-1:j+1,k]) / 9.0
    for k = 0,sz[3]-1 do sig(k) = total(sig2d[i-1:i+1,j-1:j+1,k]) / 9.0
    sig2 = stdev(sig)  ; standard deviation
    sig2 = sig2*sig2
    av = (moment(y))[0]
    chi0 = total((w*(y - av)^2)/sig2)
    chi0r = chi0/n_elements(y)
    print, 'Total chi^2 when no emission line = ', chi0
    print, 'Reduced chi^2 when no emission line = ', chi0r

    fit = curvefit(x,y,w,a,sigma,chisq=chi1,function_name="triplet",/noderivative)
    print, fit
    chi1 = total((w*(y - fit)^2)/sig2)
    chi1r = chi1/(n_elements(y)-(n_elements(a)-1))
    print, 'Total chi^2 when fit inclubed = ', chi1
    print, 'Reduced chi^2 when fit included = ', chi1r

    oplot,x,fit,color=3
    ;if (chi0-chi1 gt chi2lim) then begin

  if (chi0-chi1 gt chi2lim) then begin

     fitb = curvefit(x,y,w,b,sigma,chisq=chi1,function_name="tripletbr",/noderivative)
     chi2 = total((w*(y - fitb)^2)/sig2)
     chi2r = chi1/(n_elements(y)-(n_elements(b)-1))
     print, 'Total chi^2 when broad line inclubed = ', chi2
     print, 'Reduced chi^2 when broad line included = ', chi2r

     lim = chi2 - chi1
     if lim ge 64 then begin
       plot, wl, y
       oplot, wl, fitb, color=2
       fitcube[i,j,*] = fitb
       acube[i,j,*] = b[0:4]
       bcube[i,j,*] = b[5:7]
       print, 'best parameters = ', b
       print, '  with ', chi0-chi2, ' at ', i, j

      ;err1=errsinglet(wl,y,b,w,chi2,sig2,[0.00001,0,0,0,0,0,0,0],0,brtriplet=1)
      ;err2=errsinglet(wl,y,b,w,chi2,sig2,[0,0.001,0,0,0,0,0,0],1,brtriplet=1)
     ;err3=errsinglet(wl,y,b,w,chi2,sig2,[0,0,0.00001,0,0,0,0,0],2,brtriplet=1)
     ;err4=errsinglet(wl,y,b,w,chi2,sig2,[0,0,0,0.00001,0,0,0,0],3,brtriplet=1)
     ;err5=errsinglet(wl,y,b,w,chi2,sig2,[0,0,0,0,0.00001,0,0,0],4,brtriplet=1)
     ;err6=errsinglet(wl,y,b,w,chi2,sig2,[0,0,0,0,0,0.001,0,0],5,brtriplet=1)
     ;err7=errsinglet(wl,y,b,w,chi2,sig2,[0,0,0,0,0,0,0.0001,0],6,brtriplet=1)
     ;err8=errsinglet(wl,y,b,w,chi2,sig2,[0,0,0,0,0,0,0,0.00001],7,brtriplet=1)

     ;aerrorcube[i,j,0] = err1
     ;aerrorcube[i,j,1] = err2
     ;aerrorcube[i,j,2] = err3
     ;aerrorcube[i,j,3] = err4
     ;aerrorcube[i,j,4] = err5
     ;aerrorcube[i,j,5] = err6
     ;aerrorcube[i,j,6] = err7
     ;aerrorcube[i,j,7] = err8

     endif else begin
     plot, wl, y
     oplot, wl, fit, color=2
     fitcube[i,j,*] = fit
     acube[i,j,*] = a
     print, 'best parameters = ', a
     print, '  with ', chi0-chi1, ' at ', i, j

      ;err1=errsinglet(wl,y,a,w,chi1,sig2,[0.00001,0,0,0,0],0,triplet=1)
      ;err2=errsinglet(wl,y,a,w,chi1,sig2,[0,0.001,0,0,0],1,triplet=1)
     ;err3=errsinglet(wl,y,a,w,chi1,sig2,[0,0,0.00001,0,0],2,triplet=1)
     ;err4=errsinglet(wl,y,a,w,chi1,sig2,[0,0,0,0.00001,0],3,triplet=1)
     ;err5=errsinglet(wl,y,a,w,chi1,sig2,[0,0,0,0,0.00001],4,triplet=1)


     ;aerrorcube[i,j,0] = err1
     ;aerrorcube[i,j,1] = err2
     ;aerrorcube[i,j,2] = err3
     ;aerrorcube[i,j,3] = err4
     ;aerrorcube[i,j,4] = err5
     ;aerrorcube[i,j,5] = 0
     ;aerrorcube[i,j,6] = 0
     ;aerrorcube[i,j,7] = 0

     endelse

;wait,5

     endif else begin
        plot, wl, y
        print, ' no line detected at ', i, j
    endelse
endelse

endfor

save,finalcube,wl,acube,bcube;, aerrorcube
;writefits,'aerrorcube.fits',aerrorcube

writefits,'acube.fits',acube
writefits,'bcube.fits',bcube

stop
end

;----------------------------------------------------------------

;--------------------------------------------------------------------------------------------------------------------

pro tripletbr, x, par, f, pder

z = par[0]
zbr = par[7]
ha = 6562.8*(1+z)
NII_r = 6583*(1+z)
NII_b = 6548.1*(1+z)
inten = par[1]
ww = par[2]
contin = par[3]
NIIinten = par[4]
brint = par[5]
brww = par[6]
br = 6562.8*(1+zbr)

ww = sqrt(0.6^2 + ww^2)         ; default width to match sky line.
brww = sqrt(0.6^2 + brww^2)


;if ww ge 10 then ww = 10

Ha = (inten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-Ha)^2/ww^2))
NII_a = (NIIinten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-NII_r)^2/ww^2))
NII_b = (NIIinten/3./ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-NII_b)^2/ww^2))
bl = (brint/brww/sqrt(2.0*!pi)) * (exp(-0.5*(x-br)^2/brww^2))

f = contin + Ha + NII_a + NII_b + bl

pder = fltarr(n_elements(x),n_elements(par))   ; no value returned.

end

;----------------------------------------------------------------


pro triplet, x, par, f, pder

z = par[0]
ha = 6562.8*(1+z)
NII_r = 6583*(1+z)
NII_b = 6548.1*(1+z)
inten = par[1]
ww = par[2]
contin = par[3]
NIIinten = par[4]

ww = sqrt(0.6^2 + ww^2)         ; default width to match sky line.

Ha = (inten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-Ha)^2/ww^2))
NII_a = (NIIinten/ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-NII_r)^2/ww^2))
NII_b = (NIIinten/3./ww/sqrt(2.0*!pi)) * (exp(-0.5*(x-NII_b)^2/ww^2))

f = contin + Ha + NII_a + NII_b

pder = fltarr(n_elements(x),n_elements(par))   ; no value returned.

end
