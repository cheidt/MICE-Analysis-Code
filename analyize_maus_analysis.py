import ROOT
from config_maus_analysis import analysis_config as _config
#########################################################################################
  # 
class Analyize(object):
  def __init__ (self):
    print "INITALIZING ANALYSIS"
#    self.spaces = {"upstream":{1:{}, 2:{}, 3:{}, 4:{}, 5:{}}, \
#                 "downstream":{1:{}, 2:{}, 3:{}, 4:{}, 5:{}}}
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
          self.Draw_Line(temp)
    
  def Draw_Line(self, seeds):
    for station in seeds:
      pass
      
      
      