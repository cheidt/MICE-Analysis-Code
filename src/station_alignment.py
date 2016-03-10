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
    for detector in tracks:
      if not len(tracks[detector]) == 1:
        if len(tracks[detector]) == 0:
          _output.Message("No track in the ", detector, " tracker", exc="testing")
          continue
        else:
          _output.Message("Number of tracks in the ", detector, " ", \
                           len(tracks[detector]), exc="testing")
          pop_list = []
          for track in range(len(tracks[detector])):
            _output.Message("Track ", track, " has ", tracks[detector]\
                            [track]["triplets"], " data points", exc="testing")
            if int(tracks[detector][track]["triplets"]) <= \
               int(_config["min_trip"]):
              pop_list.append(track)
          if not len(pop_list) == 0:
            pop_list.sort(reverse=True)
            for i in range(len(pop_list)):
              tracks[detector].pop(pop_list[i])

      if not len(tracks[detector]) == 1:
        if len(tracks[detector]) == 0:
          _output.Message("No acceptable tracks", exc="testing")
          continue
        if len(tracks[detector]) > 1:
          _output.Message("Too many tracks", exc="testing")
          continue

      if tracks[detector][0]["triplets"] < _config["req_trip"]:
        _output.Message("Not enough data points in track ", \
                         tracks[detector][0]["triplets"])
        continue

      for station in tracks[detector][0]["seeds"][0][detector]:
        if not len(tracks[detector][0]["seeds"][0][detector][station]) == 1:
          _output.Message("Phantom space point in ",detector,"",station, exc="testing")

      self.spaces[detector][self.count] = {1:{}, 2:{}, 3:{}, 4:{}, 5:{}}
      for seeds in tracks[detector][0]["seeds"]:
        for station in seeds[detector]:
          temp = {}
          try:
            temp = {"x_pos":seeds[detector][station][0]["x_glob_pos"], \
                    "y_pos":seeds[detector][station][0]["y_glob_pos"], \
                    "z_pos":seeds[detector][station][0]["z_glob_pos"]}
          except IndexError:
            print seeds[detector][station]
          self.spaces[detector][self.count][station] = temp

  def Station_Alignment(self):
    for detector in self.spaces:
      for event in self.spaces[detector]:
        for cut_station in range(1,6):
          #print "event: ",event," station ",cut_station
          space = self.spaces[detector][event][cut_station]
          temp = copy.deepcopy(self.spaces[detector][event])
          temp.pop(cut_station, None)
          #print temp
          line = self.Draw_Line(temp)
          
          x_exp = line["mx"]*space["z_pos"]+line["bx"]
          x_res = x_exp - space["x_pos"]
          y_exp = line["my"]*space["z_pos"]+line["by"]
          y_res = y_exp - space["y_pos"]
          
          #print "\n",x_res,"  X Residual"
          #print y_res,"  Y Residual\n"

          name  = "st_st_residual"
          title = "Four Station Residual"
          i = "TKU" if detector == "upstream" else "TKD"
          self.o.Fill(name, title, x_res, y_res, \
                      500, -250 , 250, 500, -250 , 250, \
                      detector=i, station=cut_station)

          self.Find_Coefficents(space, x_exp, y_exp)

    self.o.Write()

  def Draw_Line(self, seeds):
    temp = {}
#    x = [1,2,3,4]
#    y = [6,5,7,10]
#    c1 = c2 = a1 = a2 = a3 = a4 = 0
#    for i in range(len(x)):
#      c1 += y[i]
#      c2 += x[i]*y[i]
#      a1 += 1
#      a2 += x[i]
#      a3 += x[i]
#      a4 += x[i]**2
#    a = np.array([[a1,a3],[a2,a4]])
#    b = np.array([c1,c2])
#    x = np.linalg.solve(a,b)
#    print x
    for dim in ["x", "y"]:
      pos = dim+"_pos"
      a1 = a2 = a3 = a4 = c1 = c2 = 0
      for station in seeds:
        c1 += seeds[station][pos]
        c2 += seeds[station][pos]*seeds[station]["z_pos"]
        a4 += seeds[station]["z_pos"]
        a3 += 1
        a2 += seeds[station][pos]**2
        a1 += seeds[station][pos]
      a = np.array([[a1,a3],[a2,a4]])
      b = np.array([c1,c2])
      #print "\n",pos
      #print "A matrix "
      #print a
      #print "B matrix"
      #print b
      x = np.linalg.solve(a,b)
      #print "X matrix"
      #print x
      m_dim = "m"+dim
      b_dim = "b"+dim
      temp[m_dim] = x[0]
      temp[b_dim] = x[1]
    return temp

  def Find_Coefficents(self, space, x_exp, y_exp):
    pass
    