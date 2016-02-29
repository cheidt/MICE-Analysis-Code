import output_maus_analysis as _output
from config_maus_analysis import analysis_config as _config
#########################################################################################
  # 
class ST_Alignment(object):
  def __init__ (self):
    print "INITALIZING ALIGNMENT"
    self.o = _output.Output("ST_Alignment")
    self.spaces = {"upstream":{}, \
                 "downstream":{}}
    self.count  = 0

  def StS_Collect_Space_Points(self, track):
    self.count += 1
    for detector in track:
      if not len(track[detector]) == 1:
        continue
      if track[detector][0]["triplets"] == 5:
        self.spaces[detector][self.count] = {1:{}, 2:{}, 3:{}, 4:{}, 5:{}}
        for seeds in track[detector][0]["seeds"]:
          for station in seeds[detector]:
            temp = {}
            temp = {"x_pos":seeds[detector][station][0]["x_glob_pos"], \
                    "y_pos":seeds[detector][station][0]["y_glob_pos"], \
                    "z_pos":seeds[detector][station][0]["z_glob_pos"]}
            self.spaces[detector][self.count][station] = temp

  def Station_Alignment(self):
    for detector in self.spaces:
      for event in self.spaces[detector]:
        for cut_station in range(1,6):
          temp = self.spaces[detector][event]
          temp.pop(cut_station, None)
          line = self.Draw_Line(temp)

  def Draw_Line(self, seeds):
    temp = {}
    for dim in ["x_pos", "y_pos"]:
      for station in seeds:
        a1 +=  2*(seeds[station]["z_pos"]**2)
        a2 += -2*(seeds[station][dim]*seeds[station]["z_pos"])
        a3 +=  2*(seeds[station]["z_pos"])
        a4 +=  8
        a5 += -2*(seeds[station][dim])
        a6 =   a3
      m_dim = "m_"+dim
      b_dim = "b_"+dim
      temp[m_dim] = (a3*a5/a4 - a2)/(a1 - a3*a6/a4)
      temp[b_dim] = -m*a1/a3 - a2/a3
    return temp