from config_maus_analysis import analysis_config as _config
import output_maus_analysis as _output
import numpy as np
import copy
#########################################################################################
  # 
class ST_Alignment(object):
  def __init__ (self):
    print "INITALIZING ALIGNMENT"
    self.o = _output.Output("station_alignment")
    self.spaces = {"upstream":{}, \
                   "downstream":{}}
    self.count  = 0

  def StS_Collect_Space_Points(self, tracks):
    self.count += 1
    print "Alignment Event: ",self.count
    for detector in tracks:
      if not len(tracks[detector]) == 1:
        if len(tracks[detector]) == 0:
          continue
        else:
          _output.Message("Number of tracks in the ", detector, " ", \
                           len(tracks[detector]))
          pop_list = []
          for track in range(len(tracks[detector])):
            _output.Message("Track ", track, " has ", tracks[detector]\
                            [track]["triplets"], " data points")
            if int(tracks[detector][track]["triplets"]) <= \
               int(_config["min_trip"]):
              pop_list.append(track)
          if not len(pop_list) == 0:
            pop_list.sort(reverse=True)
            for i in range(len(pop_list)):
              tracks[detector].pop(pop_list[i])

      if not len(tracks[detector]) == 1:
        if len(tracks[detector]) == 0:
          _output.Message("No acceptable tracks")
          continue
        if len(tracks[detector]) > 1:
          _output.Message("Too many tracks")
          continue

      if tracks[detector][0]["triplets"] < _config["req_trip"]:
        continue

      for station in tracks[detector][0]["seeds"][0][detector]:
        if not len(tracks[detector][0]["seeds"][0][detector][station]) == 1:
          _output.Message("Phantom space point in ",detector," ",station, \
                                                    exc="Always")

      for seeds in tracks[detector][0]["seeds"]:
        temp = {}
        use_event = True
        for station in seeds[detector]:
          try:
            temp[station] = {"x_pos":seeds[detector][station][0]["x_glob_pos"], \
                             "y_pos":seeds[detector][station][0]["y_glob_pos"], \
                             "z_pos":seeds[detector][station][0]["z_glob_pos"]}
          except IndexError:
            use_event = False
            continue
        if use_event == True:
          self.spaces[detector][self.count] = temp

  def Station_Alignment(self):
#    print "alignment"
    for detector in self.spaces:
      for cut_station in range(1,6):
        x_exp, y_exp, x_exp = {}, {}, {}
        for event in self.spaces[detector]:
          space = self.spaces[detector][event][cut_station]
          temp = copy.deepcopy(self.spaces[detector][event])
          temp.pop(cut_station, None)
          line = self.Draw_Line(temp)
          x_exp[event] = float(line["mx"]*space["z_pos"]+line["bx"])
          x_res        = x_exp[event] - space["x_pos"]
          y_exp[event] = float(line["my"]*space["z_pos"]+line["by"])
          y_res        = y_exp[event] - space["y_pos"]

          name  = "st_st_residual"
          title = "Four Station Residual"
          i = "TKU" if detector == "upstream" else "TKD"
          self.o.Fill(name, title, x_res, y_res, \
                      500, -250 , 250, 500, -250 , 250, \
                      detector=i, station=cut_station)

        self.Find_Coefficents(temp, x_exp, y_exp)

    self.o.Write()

  def Draw_Line(self, seeds):
    temp = {}
    for dim in ["x", "y"]:
      pos = dim+"_pos"
      a1, a2, a3, a4, c1, c2 = 0, 0, 0, 0 ,0 ,0
      for station in seeds:
        try:
          c1 += seeds[station][pos]*seeds[station]["z_pos"]
          c2 += seeds[station][pos]
          a1 += seeds[station]["z_pos"]**2
          a2 += seeds[station]["z_pos"]
          a3 += seeds[station]["z_pos"]
          a4 += 1
        except KeyError:
          print "Fit error"
          print "Station ", station
          print seeds

      a = np.matrix([[a1,a2],[a3,a4]])
      b = np.matrix([[c1],[c2]])
      x = np.linalg.solve(a,b)
      m_dim = "m"+dim
      b_dim = "b"+dim
      temp[m_dim] = x[0]
      temp[b_dim] = x[1]
    return temp

  def Find_Coefficents(self, spaces, x_exp, y_exp):
    for event in spaces:
      pass