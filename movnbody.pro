
function center,pos
;help,pos
; rarr= [particle][x, y, z]
  num=n_elements(pos)/3.0
  ret=[0.0d0,0.0d0,0.0d0]
  for i=0,num-1 do begin
    for j=0,2 do begin
      ret[j] += pos[i,j]
    endfor
  endfor
  ret /= num
  return,ret
end

pro readdata,fn,r,v,m
  openr,unit,fn,/get_lun
  num=0L & tstep=0L
  readu,unit,num,tstep
  print,num,tstep
  r=dblarr(num,3) & v=dblarr(num,3) & m=dblarr(num)

  for j=0,num-1 do begin
    rtmp=dblarr(3) & vtmp=dblarr(3) & mtmp=0.d0
    readu,unit,rtmp,vtmp,mtmp
    r[j,*] = rtmp
    v[j,*] = vtmp
    m[j] = mtmp
  endfor
  close,unit
  free_lun,unit
end


pro plotter,tstep,rad,v,m
  l=30
  cnt=center(reform(rad[*,*]))
  plot,[0,1],[0,1],/nodata ,/xstyle,/ystyle,title='t = '+numtostr(tstep) $
    ;,xrange=[0,l],yrange=[0,l]
    ;,xrange=[-l,l],yrange=[-l,l]
    ,xrange=[-l+cnt[0],l+cnt[0]],yrange=[-l+cnt[1],l+cnt[1]]
    ;,xrange=[-l+rad[0,0],l+rad[0,0]],yrange=[-l+rad[0,1],l+rad[0,1]]

  for j=0,n_elements(rad[*,0])-1 do begin
  ; oplot,rad[0:1,j,0],rad[0:1,j,1],/psym
  ; oplot,rad[0:1,j,0],rad[0:1,j,1]
    plots,circle(rad[j,0],rad[j,1],0.1)
  endfor

  ;wait,0.0001
  pngfn='../img/img'+numtostr(tstep,format='(I09)')+'.png'
  tvlct,r,g,b,/get
  write_png,pngfn,tvrd(),r,g,b
end


;window,xsize=600,ysize=600
set_plot,'Z'
device,set_resolution=[600,600]
rainbowf
r=0
v=0
m=0
;readdata,"./dat/st00000000.bin",rad,v,m
i=0
DataTop="../dat/"
readdata,DataTop+"/st"+numtostr(i,format='(I08)')+".bin",rad,v,m
plotter,0,rad,v,m


fnlist=file_search(DataTop+"/*bin")
for i=1,n_elements(fnlist)-1 do begin
  readdata,DataTop+"/st"+numtostr(i,format='(I08)')+".bin",rad,v,m
  plotter,i,rad,v,m

endfor

device,/close
end
