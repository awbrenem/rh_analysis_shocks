;Dombeck's program based on Lynn's rh_crib_example.txt


pro rh_analysis,bdat,vdat,pdat,tdat,ttl=ttl

  bbdat=bdat
  vvdat=vdat
  ppdat=pdat
  ttdat=tdat

  badb=where(bbdat[*,1:3] lt -1.e29,cntb)
  badv=where(vvdat[*,1:3] lt -1.e29,cntv)
  badp=where(ppdat[*,1] lt -1.e29,cntp)
  badt=where(ttdat[*,1] lt -1.e29,cntt)
  if cntb gt 0 then bbdat[badb,1:3]=!values.d_nan
  if cntv gt 0 then vvdat[badv,1:3]=!values.d_nan
  if cntp gt 0 then ppdat[badp,1]=!values.d_nan
  if cntt gt 0 then ttdat[badt,1]=!values.d_nan

  totplot,bbdat,[['UT','R-H Bx','R-H By','R-H Bz'],['sec','nT','nT','nT']]
  totplot,vvdat,[['UT','R-H Vx','R-H Vy','R-H Vz'],['sec','km/s','km/s','km/s']]
  totplot,ppdat,[['UT','R-H Density'],['sec','cm^-3']]
  totplot,ttdat,[['UT','R-H Temperature'],['sec','Energy']]

  tplot,['r-hbx','r-hby','r-hbz','r-hvx','r-hvy','r-hvz','r-hdensity','r-htemperature']

select:
  dn=0
  while dn eq 0 do begin
    print,'Enter Zoom: 1=In, 2=Out, 3=All, 0=Done'
    inp=''
    read,inp
      case inp of
       '1' : tzoom
       '2' : tzoom,/out
       '3' : tzoom,/all
       '0' : dn=1
       else : print,'Enter 0-3'
      endcase
  endwhile

  dn=0
  status=0
  while dn eq 0 do begin
    print,'Select Upstream Time Range: 1=Select, 2=Start over at Zoom. 0=Done'
    inp=''
    read,inp
      case inp of
       '1' : status=select_times(upt1,upt2)
       '2' : goto,select
       '0' : if status eq 1 then dn=1
       else : print,'Enter 0-1'
      endcase
  endwhile

  dn=0
  status=0
  while dn eq 0 do begin
    print,'Select Downstream Time Range: 1=Select, 2=Start over at Zoom, 0=Done'
    inp=''
    read,inp
      case inp of
       '1' : status=select_times(dwt1,dwt2)
       '2' : goto,select
       '0' : if status eq 1 then dn=1
       else : print,'Enter 0-1'
      endcase
  endwhile

  upbbdat=time_trim(bbdat,upt1,upt2)
  upvvdat=time_trim(vvdat,upt1,upt2)
  upppdat=time_trim(ppdat,upt1,upt2)
  upttdat=time_trim(ttdat,upt1,upt2)
  upminpts=min([n_elements(upbbdat[*,0]),n_elements(upvvdat[*,0]),n_elements(upppdat[*,0]),n_elements(upttdat[*,0])],upmnix)
  if upminpts le 1 then begin
    print,'Not enough Upstream Points'
    goto,select
  endif
    case upmnix of
      0 : begin
            upbdat=upbbdat
            upvdat=resample(upvvdat,upbbdat[*,0],/skipbad)
            uppdat=resample(upppdat,upbbdat[*,0],/skipbad)
            uptdat=resample(upttdat,upbbdat[*,0],/skipbad)
          end
      1 : begin
            upvdat=upvvdat
            upbdat=resample(upbbdat,upvvdat[*,0],/skipbad)
            uppdat=resample(upppdat,upvvdat[*,0],/skipbad)
            uptdat=resample(upttdat,upvvdat[*,0],/skipbad)
          end
      2 : begin
            uppdat=upppdat
            upvdat=resample(upvvdat,upppdat[*,0],/skipbad)
            upbdat=resample(upbbdat,upppdat[*,0],/skipbad)
            uptdat=resample(upttdat,upppdat[*,0],/skipbad)
          end
      3 : begin
            uptdat=upttdat
            upvdat=resample(upvvdat,upttdat[*,0],/skipbad)
            uppdat=resample(upppdat,upttdat[*,0],/skipbad)
            upbdat=resample(upbbdat,upttdat[*,0],/skipbad)
          end
    endcase

  dwbbdat=time_trim(bbdat,dwt1,dwt2)
  dwvvdat=time_trim(vvdat,dwt1,dwt2)
  dwppdat=time_trim(ppdat,dwt1,dwt2)
  dwttdat=time_trim(ttdat,dwt1,dwt2)
  dwminpts=min([n_elements(dwbbdat[*,0]),n_elements(dwvvdat[*,0]),n_elements(dwppdat[*,0]),n_elements(dwttdat[*,0])],dwmnix)
  if dwminpts le 1 then begin
    print,'Not enough Upstream Points'
    goto,select
  endif
    case dwmnix of
      0 : begin
            dwbdat=dwbbdat
            dwvdat=resample(dwvvdat,dwbbdat[*,0],/skipbad)
            dwpdat=resample(dwppdat,dwbbdat[*,0],/skipbad)
            dwtdat=resample(dwttdat,dwbbdat[*,0],/skipbad)
          end
      1 : begin
            dwvdat=dwvvdat
            dwbdat=resample(dwbbdat,dwvvdat[*,0],/skipbad)
            dwpdat=resample(dwppdat,dwvvdat[*,0],/skipbad)
            dwtdat=resample(dwttdat,dwvvdat[*,0],/skipbad)
          end
      2 : begin
            dwpdat=dwppdat
            dwvdat=resample(dwvvdat,dwppdat[*,0],/skipbad)
            dwbdat=resample(dwbbdat,dwppdat[*,0],/skipbad)
            dwtdat=resample(dwttdat,dwppdat[*,0],/skipbad)
          end
      3 : begin
            dwtdat=dwttdat
            dwvdat=resample(dwvvdat,dwttdat[*,0],/skipbad)
            dwpdat=resample(dwppdat,dwttdat[*,0],/skipbad)
            dwbdat=resample(dwbbdat,dwttdat[*,0],/skipbad)
          end
    endcase

    npts=min([upminpts,dwminpts])
    if upminpts lt dwminpts then begin
      dwbdat=dwbdat[0:upminpts-1,*]
      dwvdat=dwvdat[0:upminpts-1,*]
      dwpdat=dwpdat[0:upminpts-1,*]
      dwtdat=dwtdat[0:upminpts-1,*]
    endif else if upminpts gt dwminpts then begin
      upbdat=upbdat[0:dwminpts-1,*]
      upvdat=upvdat[0:dwminpts-1,*]
      uppdat=uppdat[0:dwminpts-1,*]
      uptdat=uptdat[0:dwminpts-1,*]
    endif

  tup1=time_string(upbdat[0,0],/msec)
  tup2=time_string(upbdat[npts-1,0],/msec)
  tdw1=time_string(dwbdat[0,0],/msec)
  tdw2=time_string(dwbdat[npts-1,0],/msec)
  tdate=strmid(tup1,0,10)
  tupst='["'+strmid(tup1,11,12)+'"],["'+strmid(tup2,11,12)+'"]'
  tdwst='["'+strmid(tdw1,11,12)+'"],["'+strmid(tdw2,11,12)+'"]'

  upbdat=upbdat[*,1:3]
  upvdat=upvdat[*,1:3]
  uppdat=uppdat[*,1]
  uptdat=uptdat[*,1]

  dwbdat=dwbdat[*,1:3]
  dwvdat=dwvdat[*,1:3]
  dwpdat=dwpdat[*,1]
  dwtdat=dwtdat[*,1]

;; avsw   = [[[TRANSPOSE(vswup)]],[[TRANSPOSE(vswdn)]]]        ; => [N,3,2]-Element Array
;; amagf  = [[[TRANSPOSE(magup)]],[[TRANSPOSE(magdn)]]]        ; => [N,3,2]-Element Array
;; adens  = [[denup],[dendn]]                                  ; => [N,2]-Element Array
;; atemp  = [[tempup],[tempdn]]                                ; => [N,2]-Element Array

   avsw=[[[upvdat]],[[dwvdat]]]
   amagf=[[[upbdat]],[[dwbdat]]]
   adens=[[uppdat],[dwpdat]]
   atemp=[[uptdat],[dwtdat]]

stop

; => Generate dummy array of angles
m          = 100L
phi        = DINDGEN(m)*2d0*!DPI/(m - 1L)
the        = DINDGEN(m)*!DPI/(m - 1L) - !DPI/2d0
ph         = REFORM(phi)
th         = REFORM(the)
; => Generate shock normal vector
;            [theta, phi, 3]
nor        = DBLARR(m,m,3L)
nor[*,*,0] = COS(th) # COS(ph)
nor[*,*,1] = COS(th) # SIN(ph)
nor[*,*,2] = SIN(th) # REPLICATE(1,m)
;-----------------------------------------------------------------------------------------
; => Calculate chi-squared distribution using Equations 2, 3, 4, 5, and 6 for 1997-12-10
;-----------------------------------------------------------------------------------------
nqq        = [1,1,1,1,1]
chisq      = rh_solve_lmq(adens,avsw,amagf,atemp,NEQS=nqq,SOLN=soln)
print
print,'Date: '+tdate
print,'Upstream Times: '+tupst
print,'Dwstream Times: '+tdwst
print,'# of points: ',+npts
print
; => Print out best fit angles
PRINT,'Best Fit/Std Dev - Theta: ', soln.THETA*18d1/!DPI
;         Avg.          Std. Dev.
;      -22.290909       6.3262974
PRINT,'Best Fit/Std Dev - Phi: ', soln.PHI*18d1/!DPI
;         Avg.          Std. Dev.
;       167.56364       2.9541958

; => Print out best fit shock normal speed in spacecraft frame [km/s]
PRINT,'Best Fit/Std Dev - SC-Frame Speed(km/s): ', soln.VSHN
;       401.05219       27.929917

; => Print out best fit upstream shock normal speed in shock frame [km/s]
PRINT,'Best Fit/Std Dev - SH-Frame Speed(km/s) ', soln.USHN_UP
;      -133.18804       16.392315

; => Print out best fit shock normal vector [GSE coordinates]
PRINT,'Best Fit - Normal Vector: ', soln.SH_NORM[*,0]
;     -0.89712733      0.19784031     -0.37712047

; => Print out uncertainty of shock normal vector [GSE coordinates]
PRINT,'Uncertainty - Normal Vector', soln.SH_NORM[*,1]
;     0.040895515     0.046225610      0.10007386

stop

;-----------------------------------------------------------------------------------------
; => Region of 68.3% confidence -> ∆X^2 = 2.30 (for 2 degrees of freedom)
;      Theorem D from numerical recipes [Section 15.6]
;-----------------------------------------------------------------------------------------
mnchsq  = MIN(chisq,/NAN,ln)
conf683 = mnchsq[0] + 2.30
region  = WHERE(chisq LE conf683[0],greg)
rind    = ARRAY_INDICES(chisq,region)
;chi_reg = REPLICATE(d,100L,100L)
;chi_reg[rind[0,*],rind[1,*]] = chisq[rind[0,*],rind[1,*]]
;-----------------------------------------------------------------------------------------
; => plot chi-squared
;-----------------------------------------------------------------------------------------
yra       = [-9d1,9d1]
xra       = [0d0,36d1]
xttl      = 'phi (deg)'
yttl      = 'theta (deg)'

nlevels   = 60L
result    = ALOG(chisq*1d0)/ALOG(1d1)
range     = [MIN(result,/NAN),MAX(result,/NAN)]
cpd0      = FLOOR(nlevels/(range[1] - range[0])) > 1.
cpd       = FIX(nlevels/(range[1] - range[0])) > 1
nn2       = nlevels - 1L
mxres     = CEIL(MAX(result,/NAN))
levels    = FINDGEN(nlevels)*(mxres[0] - range[0])/nn2 + range[0]
levels    = roundsig(levels,SIG=2)
color     = LONARR(nlevels)
color     = LINDGEN(nlevels)*(254L - 15L)/(nlevels - 1L) + 15L
color     = ROUND((color*1e3)^(2e0/3e0)/16e0)
c_colors  = BYTE(color)
good_lbs           = INDGEN(nlevels/3L)*3L + 2L
c_labels           = REPLICATE(0,nlevels)
c_labels[good_lbs] = 1

ttle      = '!7v!3'+'!U2!N'+'(!7h!3'+',!7u!3'+') for an interplanetary shock on: '+tdate[0]
pstr      = {NLEVELS:nlevels,XTITLE:xttl,YTITLE:yttl,$
             XSTYLE:1,YSTYLE:1,ZSTYLE:1,FILL:1,C_COLORS:c_colors,   $
             XRANGE:xra,YRANGE:yra,XLOG:0,YLOG:0,ZLOG:0,XMINOR:11L, $
             YMINOR:11L,TITLE:ttle}

WINDOW,1,RETAIN=2
WSET,1
!P.MULTI = 0
CONTOUR,TRANSPOSE(result),phi*18d1/!DPI,the*18d1/!DPI,_EXTRA=pstr
  OPLOT,[soln.PHI[0]*18d1/!DPI],[soln.THETA[0]*18d1/!DPI],PSYM=2,COLOR=250
  ; => Plot "ellipse" containing the 68.3% confidence interval
  CONTOUR,TRANSPOSE(result),phi*18d1/!DPI,the*18d1/!DPI,/OVERPLOT,XSTYLE=0,YSTYLE=0,$
          C_THICK=2,C_LINESTYLE=3,LEVELS=ALOG(conf683[0])/ALOG(1d1)
wset,0

end
