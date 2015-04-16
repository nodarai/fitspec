pro countzeroes, finalcube, value

  zeroes=0
  sz=size(finalcube)
  for i=0,sz[1]-1 do for j=0,sz[2]-1 do for k=0,sz[3]-1 do begin
    if long(finalcube[i,j,k]) eq long(value) then zeroes=zeroes+1
  endfor
  print, zeroes

end
