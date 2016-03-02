import output_maus_analysis as _output
from config_maus_analysis import analysis_config as _config

"""
This script holds low level analysis methods, mostly intented to print out
quick histograms without doing to much heavy lifting. 
"""


class Analysis(object):
  def __init__ (self):
    self.o_sp = _output.Output("space_point")
    self.o_TtT = _output.Output("TOF_to_TOF")
    print "INITALIZING ANALYSIS"

#########################################################################################
  # Compares Truth position to reconstructed positions
  def SP_to_Virt(self, virtuals, spaces):
    for detector in spaces:
      for station in spaces[detector]:
        if len(spaces[detector][station]) == 1:
          for virtual in virtuals[detector][station]:
            if virtual["plane"] == 0:
              space   = spaces[detector][station][0]
              x_res   = virtual["x_pos"] - space["x_glob_pos"]
              y_res   = virtual["y_pos"] - space["y_glob_pos"]
              i = "TKU" if detector == "upstream" else "TKD"
              name  = "SP_to_Virt"
              title = "Residuals MC Truth to Space Point" 
              self.o.Fill(name, title, x_res, y_res, \
                          500, -10, 10, 500, -10, 10, \
                          detector=i, station=station)
              name = "SP_to_Virt_Hist_x"
              title = "Residuals MC Truth to Space Point X"
              self.o.Fill(name, title, x_res, 500, -10, 10, \
                          detector=i, station=station)
              name = "SP_to_Virt_Hist_x"
              title = "Residuals MC Truth to Space Point Y"
              self.o.Fill(name, title, x_res, 500, -10, 10, \
                          detector=i, station=station)

#########################################################################################
  # Pulls out a plots tracker space point positions
  def SP_Fill_ROOT(self, space_points):
    for detector in space_points:
      for station in space_points[detector]:
        if len(space_points[detector][station]) == 1:
          x_pos   = space_points[detector][station][0]["x_glob_pos"]
          y_pos   = space_points[detector][station][0]["y_glob_pos"]
          i = "TKU" if detector == "upstream" else "TKD"
          name  = "SP_Pos"
          title = "Space Point Positions"
          self.o_sp.Fill(name, title, x_pos, y_pos, \
                         250, -250 , 250, 250, -250 , 250, \
                         detector=i, station=station)
          if space_points[detector][station][0]["type"] == 3:
            name  = "SP_Trip_Pos"
            title = "Space Point Triplet Postions"
            self.o_sp.Fill(name, title, x_pos, y_pos, \
                           250, -250 , 250, 250, -250 , 250, \
                           detector=i, station=station)
          else:
            name  = "SP_Doub_Pos"
            title = "Space Point Doublet Postions"
            self.o_sp.Fill(name, title, x_pos, y_pos, \
                           250, -250 , 250, 250, -250 , 250, \
                           detector=i, station=station)

#########################################################################################
  # Writes out MC truth positions, catorigizes virtual planes near tracker
  #   stations as virtual "tracker hits"
  def Virt_Fill_ROOT(self, virtuals):
    for detector in virtuals:
      if detector == "uncat":
        for virtual in virtuals[detector]:
          x_pos = virtual["x_pos"]
          y_pos = virtual["y_pos"]
          name  = "Virt_Pos"
          title = "Truth Postions"
          self.o.Fill(name, title, x_pos, y_pos, \
                      250, -250 , 250, 250, -250 , 250, \
                      detector=detector)
      else:
        for station in virtuals[detector]:
          for virtual in virtuals[detector][station]:
            x_pos   = virtual["x_pos"]
            y_pos   = virtual["y_pos"]
            i = "TKU" if detector == "upstream" else "TKD"
            name  = "Virt_Pos"
            title = "Truth Postions"
            self.o.Fill(name, title, x_pos, y_pos, \
                        250, -250 , 250, 250, -250 , 250, \
                        detector=i, station=station)

#########################################################################################
  # Produces timing information between two TOFS
  def TOF_Timing_Info(self, up_tof, down_tof):
    if len(up_tof) == 1 and \
       len(down_tof) == 1:
      timing = up_tof[0]["time"] - down_tof[0]["time"]
      return timing
    else:
      timing = -10000.
      return timing

#########################################################################################
  # Reads in TOF1 and TOF2 positions and draws a line between the two.  Space
  #   points are transformed into global coordinates and the residuals between
  #   where the TOF to TOF line cross the tracker plane and space points are 
  #   calculated.
  def TOF_to_TOF_Tkr_Res(self, tracks, tof1, tof2):
    for detector in tracks:
      if tracks[detector]["triples"] < _config["TtT_trip_req"]:
        continue
      TOF_distance = tof1["z_pos"] - tof2["z_pos"]
      x_change = tof1["x_pos"] - tof2["x_pos"]
      y_change = tof1["y_pos"] - tof2["y_pos"]
      x_slope  = x_change/TOF_distance
      y_slope  = y_change/TOF_distance

      for spaces in track[detector]["seeds"]:
        for station in spaces[detector]:
          space = spaces[detector][station]
          z_pos = space["z_glob_pos"]
          distance = z_pos - tof1["z_pos"]
          expected_x = tof1["x_pos"] + distance * x_slope
          residual_x = space["x_glob_pos"] - expected_x
          expected_y = tof1["y_pos"] + distance * y_slope
          residual_y = space["y_glob_pos"] - expected_y
          
          name  = "TOF_to_TOF"
          title = "TOF Line and Track SP Residuals"
          self.o_TtT.Fill(name, title, x_pos, y_pos, \
                          500, -400 , 400, 500, -400 , 400, \
                          detector=i, station=station)
          
          name  = "TOF_to_TOF_X"
          title = "TOF Line and Track SP Residuals X"
          self.o_TtT.Fill(name, title, x_pos, -400 , 400, \
                          detector=i, station=station)
          
          name  = "TOF_to_TOF_Y"
          title = "TOF Line and Track SP Residuals Y"
          self.o_TtT.Fill(name, title, y_pos, 500, -400 , 400, \
                          detector=i, station=station)

  def Write(self):
    self.o_TtT.Write()
    self.o_sp.Write()