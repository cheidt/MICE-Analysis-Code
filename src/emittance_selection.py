import math

def Diff_Scrap(track_pt):   
  x_i = track_pt["x_pos"]
  y_i = track_pt["y_pos"]
  x_p = track_pt["x_mom"]
  y_p = track_pt["y_mom"]
  z_p = track_pt["z_mom"]
  d_z = 800 #mm

  x_f = (d_z*x_p/z_p)+x_i
  y_f = (d_z*y_p/z_p)+y_i
  radius = math.sqrt(x_f**2 + y_f**2)

  if radius > 90:
    return True
  return False