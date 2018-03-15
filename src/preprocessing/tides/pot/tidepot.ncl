begin
  Lp = 760
  Mp = 630
;  Lp = 542
;  Mp = 770

  pots = asciiread("waom10_8con.dat", (/Mp,Lp,16/), "float")
;  pots = asciiread("Bering_8con.dat", (/Mp,Lp,16/), "float")


  fname = "waom10_ptds.nc"
;  fname = "/archive/u1/uaf/kate/Arctic/Files/Arctic_tides_otps.nc"
  ncid = addfile(fname, "c")
  Ntc = 8
  tide_Pamp = new((/Ntc, Mp, Lp/), double)
  tide_Pphase = new((/Ntc, Mp, Lp/), double)

; Fill in the values
; ROMS order is q1, o1, p1, k1, n2, m2, s2, k2 while the input file is 
; q1, o1, p1, k1, n2, m2, s2, k2 in amplitude, phase pairs
  tide_Pamp(0,:,:) = pots(:,:,0)
  tide_Pamp(1,:,:) = pots(:,:,2)
  tide_Pamp(2,:,:) = pots(:,:,4)
  tide_Pamp(3,:,:) = pots(:,:,6)
  tide_Pamp(4,:,:) = pots(:,:,8)
  tide_Pamp(5,:,:) = pots(:,:,10)
  tide_Pamp(6,:,:) = pots(:,:,12)
  tide_Pamp(7,:,:) = pots(:,:,14)
  tide_Pphase(0,:,:) = pots(:,:,1)
  tide_Pphase(1,:,:) = pots(:,:,3)
  tide_Pphase(2,:,:) = pots(:,:,5)
  tide_Pphase(3,:,:) = pots(:,:,7)
  tide_Pphase(4,:,:) = pots(:,:,9)
  tide_Pphase(5,:,:) = pots(:,:,11)
  tide_Pphase(6,:,:) = pots(:,:,13)
  tide_Pphase(7,:,:) = pots(:,:,15)

  tide_Pamp!0 = "tide_period"
  tide_Pamp!1 = "eta_rho"
  tide_Pamp!2 = "xi_rho"
  tide_Pphase!0 = "tide_period"
  tide_Pphase!1 = "eta_rho"
  tide_Pphase!2 = "xi_rho"
  
  tide_Pamp@units = "meter"
  tide_Pamp@long_name = "tidal potential elevation amplitude"
  tide_Pamp@field = "POT_Tamp, scalar, series"
  tide_Pphase@units = "degrees, time of maximum elevation with respect chosen time origin"
  tide_Pphase@long_name = "tidal potential elevation phase angle"
  tide_Pphase@field = "POT_Tphase, scalar, series"

; global attributes
  fileAtt = True
  fileAtt@history = "adding potential tides"
  fileattdef(ncid, fileAtt)

; Write to file
  ncid->tide_Pamp = tide_Pamp
  ncid->tide_Pphase = tide_Pphase

end
