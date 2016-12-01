import output_maus_analysis as _output
import analyize_maus_analysis as _analysis
import emittance_selection as _selection
from config_maus_analysis import emittance_config as _config
import math
import numpy as np

class Analysis(object):
  def __init__ (self):
    self.o_sel = _output.Output("Selection")
    _output.Message("INITALIZING EMITTANCE", exc=True)
    self.terms = {"upstream":[], \
                "downstream":[]}
    self.event_cut = {"upstream":{"Diff_Scrap":[]}, \
                    "downstream":{"Diff_Scrap":[]}}

  def Process (self, data, timing):
    _output.Message("Collecting Emittance Data")
#    self.Plot()
    for detector in data:
      tracks = data[detector]
      for i in range(len(tracks)):
        for j in range(len(tracks[i]["track_points"])):
          track_pt = tracks[i]["track_points"][j]
          if track_pt["station"] == 5 and track_pt["plane"] == 0:
            _output.Message("Found a {} point".format(detector))
            self.terms[detector].append([track_pt["x_pos"], \
                                         track_pt["y_pos"], \
                                         track_pt["z_pos"], \
                                         track_pt["x_mom"], \
                                         track_pt["y_mom"], \
                                         track_pt["z_mom"]])
            if _config["cut_list"]["Diff_Scrap"] == True:
              self.event_cut[detector]["Diff_Scrap"].\
                             append(_selection.Diff_Scrap(track_pt))

#  def Plot(self):
#    _output.print_data(self.data)
#    for detector in self.data:
#      tracks = self.data[detector]
#      for i in range(len(tracks)):
#        for j in range(len(tracks[i]["track_points"])):
#          track_pt = tracks[i]["track_points"][j]

#          if track_pt["station"] == 1 and track_pt["plane"] == 0:
#            mom = track_pt["z_mom"]
#            tof_time = self.timing["upstream"]
#            name = "MvT"
#            title = "Momentum v Time"
#            self.o_sel.Fill(name, title, mom, tof_time, \
#                            250, 100, 300, 250, 26, 36, \
#                            detector = detector)
              
#          if track_pt["station"] == 1 and track_pt["plane"] == 0 and \
#             track_pt["scrapping"] == True:
#            mom = track_pt["z_mom"]
#            tof_time = self.timing["upstream"]
#            name = "MvTwDifCut"
#            title = "Momentum v Time with Disffuser Cut"
#            self.o_sel.Fill(name, title, mom, tof_time, \
#                            250, 100, 300, 250, 26, 36, \
#                            detector = detector)

  def Write(self):
    self.o_sel.Write()
    
  def Emittance(self):
    for detector in self.terms:
      data = np.array(self.terms[detector]).T
      cov  = np.cov(data)
      print detector, '\n ', cov, '\n'
      emit = abs(float(np.linalg.det(cov)))**(1./6.)
      print emit, '\n'