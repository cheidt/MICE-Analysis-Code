import output_maus_analysis as _output
from config_maus_analysis import analysis_config as _config
#########################################################################################
  # 
class Analysis(object):
  def __init__ (self):
    print "INITALIZING ANALYSIS"
    self.o =  _output.Output("Analysis")

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
                          
  def Write(self):
    self.o.Write()