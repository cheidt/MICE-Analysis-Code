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
      timing = down_tof[0]["time"] - up_tof[0]["time"]
      return timing
    else:
      timing = -10000.
      return timing

#########################################################################################
  # This method does some additional data cutting and then calls different
  #   methods to determine the residuals of TOF to Tracker reconstruction. 
  def TOF_Tk_Res(self, tracks, tof1, tof2):
    for detector in tracks:
      if not len(tracks[detector]) == 1:
        if len(tracks[detector]) == 0:
          continue
        else:
          print "Number of tracks in Tk ",detector," ",len(tracks[detector])
          pop_list = []
          for track in range(len(tracks[detector])):
            print "Track ",track," has ",tracks[detector][track]\
                  ["number_data_points"]," data points"
            if int(tracks[detector][track]["number_data_points"]) <= \
               int(_config["min_data"]):
              pop_list.append(track)
          pop_list.sort(reverse=True)
          if not len(pop_list) == 0:
            for i in range(len(pop_list)):
              tracks[detector].pop(pop_list[i])
          if not len(tracks[detector]) == 1:
            if len(tracks[detector]) == 0:
              print "No acceptable tracks"
              continue
            if len(tracks[detector]) > 1:
              print "Still too many tracks"
              continue
      if tracks[detector][0]["number_data_points"] < _config["tk_data_req"]:
        continue
      if tracks[detector][0]["p_value"] < 0.01:
        continue
      if not len(tracks[detector][0]["prtrks"]) == 1:
        print "Wrong number of PR tracks in track"
        continue
      if not len(tof1) == 1 or not len(tof2) == 1:
        print "Event has too many TOF spacepoints: ","\nTOF1: ",len(tof1), \
              "\nTOF2: ", len(tof2)
        continue
      tof = tof1 if detector == "upstream" else tof2
      self.TOF_to_TOF_Tk_Res(tracks[detector][0], tof1, tof2, detector)
      self.Tracker_to_TOF_Res(tracks[detector][0]["prtrks"][0], tof, detector)

#########################################################################################
  # Draws a line between TOF1 and TOF2 then calculates where in x/y that line
  #   crosses each station.  Residual is between that point and reconstructed
  #   station space point position.
  def TOF_to_TOF_Tk_Res(self, track, tof1, tof2, detector):
      TOF_distance = tof1[0]["z_pos"] - tof2[0]["z_pos"]
      x_change = tof1[0]["x_pos"] - tof2[0]["x_pos"]
      y_change = tof1[0]["y_pos"] - tof2[0]["y_pos"]
      x_slope  = x_change/TOF_distance
      y_slope  = y_change/TOF_distance

      for point in track["track_points"]:
        if point["plane"] == 0:
          z_pos = point["z_pos"]
          distance = z_pos - tof1[0]["z_pos"]
          expected_x = tof1[0]["x_pos"] + distance * x_slope
          residual_x = point["x_pos"] - expected_x
          expected_y = tof1[0]["y_pos"] + distance * y_slope
          residual_y = point["y_pos"] - expected_y

          station = point["station"]
          name  = "TOF_to_TOF"
          title = "TOF Line to Trk SP Residuals"
          self.o_TtT.Fill(name, title, residual_x, residual_y, \
                          500, -400 , 400, 500, -400 , 400, \
                          detector=detector, station=station)

          name  = "TOF_to_TOF_X"
          title = "TOF Line to Trk SP Residuals X"
          self.o_TtT.Fill(name, title, residual_x, 500, -400 , 400, \
                          detector=detector, station=station)

          name  = "TOF_to_TOF_Y"
          title = "TOF Line to Trk SP Residuals Y"
          self.o_TtT.Fill(name, title, residual_y, 500, -400 , 400, \
                          detector=detector, station=station)

#########################################################################################
  # Uses straight PR line and traces where line crosses either TOF1 or TOF2
  #   x/y values at the point are found and residuals with TOF space point 
  #   calculated.
  def Tracker_to_TOF_Res(self, track, tof, detector):
    exp_x      = track["m_x"]*tof[0]["z_pos"]+track["x_0"]
    exp_y      = track["m_y"]*tof[0]["z_pos"]+track["y_0"]
    residual_x = exp_x - tof[0]["x_pos"]
    residual_y = exp_y - tof[0]["y_pos"]

    name  = "Tk_to_TOF"
    title = "Tk Line to TOF SP Residuals"
    self.o_TtT.Fill(name, title, residual_x, residual_y, \
                    500, -400 , 400, 500, -400 , 400, \
                    detector=detector)

    name  = "Tk_to_TOF_X"
    title = "Tk Line to TOF SP Residuals X"
    self.o_TtT.Fill(name, title, residual_x, 500, -400 , 400, \
                    detector=detector)

    name  = "Tk_to_TOF_Y"
    title = "Tk Line to TOF SP Residuals Y"
    self.o_TtT.Fill(name, title, residual_y, 500, -400 , 400, \
                    detector=detector)

#########################################################################################
  # Writes analysis out to file by calling output class
  def Write(self):
    self.o_TtT.Write()
    self.o_sp.Write()