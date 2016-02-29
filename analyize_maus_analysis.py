import output_maus_analysis as _output
from config_maus_analysis import analysis_config as _config

"""
This script holds low level analysis methods, mostly intented to print out
quick histograms without doing to much heavy lifting. 
"""


class Analysis(object):
  def __init__ (self):
    print "INITALIZING ANALYSIS"
    self.o = _output.Output("analysis")

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
          self.o.Fill(name, title, x_pos, y_pos, \
                      250, -250 , 250, 250, -250 , 250, \
                      detector=i, station=station)
          if space_points[detector][station][0]["type"] == 3:
            name  = "SP_Trip_Pos"
            title = "Space Point Triplet Postions"
            self.o.Fill(name, title, x_pos, y_pos, \
                        250, -250 , 250, 250, -250 , 250, \
                        detector=i, station=station)
          else:
            name  = "SP_Doub_Pos"
            title = "Space Point Doublet Postions"
            self.o.Fill(name, title, x_pos, y_pos, \
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
            self.o.Fill(name, title, x_pos, y_pos, \
                        250, -250 , 250, 250, -250 , 250, \
                        detector=i, station=station)

#########################################################################################
  # Produces timing information between two TOFS
  def TOF_Timing_Info(self, up_tof, down_tof):
    if len(up_tof) == 1 and \
       len(down_tof) == 1:
      timing = up_tof["time"] - down_tof["time"]
      return timing
    else:
      timing = -10000.
      return timing


  def Write(self):
    self.o.Write()