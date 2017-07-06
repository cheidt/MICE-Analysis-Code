import output_maus_analysis as _output
import numpy as np
import json

class Analysis(object):
  def __init__ (self):
    self.o_sel = _output.Output("Selection")
    _output.Message("INITALIZING EMITTANCE", exc=True)
    self.terms = []
    self.data = []
    self.factor = {"upstream":[], "downstream":[]}
    data = open('emitt_data', 'w')
    data.close()

  def Process (self, data, cut_list):
    _output.Message("Collecting Emittance Data")
    event = {"upstream":[], "downstream":[]}
    for detector in data:
      tracks = data[detector]
      for i in range(len(tracks)):
        for j in range(len(tracks[i]["track_points"])):
          track_pt = tracks[i]["track_points"][j]
          if track_pt["station"] == 5 and track_pt["plane"] == 0 and \
                                  track_pt["x_mom"]**2 > .000001  and track_pt["y_mom"]**2 > .000001:
            _output.Message("Found a {} point".format(detector))
            if cut_list["Pass"]:
              event[detector].append([track_pt["x_pos"], \
                                      track_pt["x_mom"], \
                                      track_pt["y_pos"], \
                                      track_pt["y_mom"]])
            #                                    track_pt["z_mom"]])
            #            event[detector].append([track_pt["x_pos"], \
            #                                    track_pt["y_pos"], \
            #                                    track_pt["z_pos"], \
            #                                    track_pt["x_mom"], \
            #                                    track_pt["y_mom"], \
              self.factor[detector].append(track_pt["z_mom"]/105.66)
              name = "Emit_X_{}".format(detector)
              title = "Emittance X {}".format(detector)
              self.o_sel.Fill(name, title, track_pt["x_pos"], track_pt["x_mom"], \
                              300, -200, 200, 300, -200, 200)
              name = "Emit_Y_{}".format(detector)
              title = "Emittance Y {}".format(detector)
              self.o_sel.Fill(name, title, track_pt["y_pos"], track_pt["y_mom"], \
                              300, -200, 200, 300, -200, 200)
              self.terms.append(event)
              #self.terms.append([event, cut_list])

  def Write(self):
    self.o_sel.Write()

  #########################################################################################
  # Takes in emittance terms and outputs system emittance.  Separates out cut booleans
  # and uses those to determine what makes it to final calculation.  Solves for
  # 4D and the two 2D emittances.
  def Emittance(self):
    eval   = {"upstream": [], "downstream": []}
    x_eval = {"upstream": [], "downstream": []}
    y_eval = {"upstream": [], "downstream": []}
    passed = {"upstream": 0, "downstream": 0}
    count  = {"upstream": 0, "downstream": 0}
    for event in self.terms:
#      cut_list = event[-1]
#      del event[-1]
      for detector in event:
        for i in range(len(event[detector])):
          x_array = []
          y_array = []
          eval[detector].append(event[detector][i])
          for dim in range(len(event[detector][i])):
            if dim < 2:
              x_array.append(event[detector][i][dim])
            else:
              y_array.append(event[detector][i][dim])
          x_eval[detector].append(x_array)
          y_eval[detector].append(y_array)
    for detector in eval:
      data   = np.array(eval[detector]).T
      x_data = np.array(x_eval[detector]).T
      y_data = np.array(y_eval[detector]).T
      #      for dim in range(len(data)):
      #        ave = sum(data[dim])/len(data[dim])
      #        for num in range(len(data[dim])):
      #          data[dim][num] = data[dim][num] - ave
      cov    = np.cov(data)
      x_cov  = np.cov(x_data)
      y_cov  = np.cov(y_data)
      print detector
      print "\nX Y Covariance:"
      print cov, "\n"
      emit   = sum(self.factor[detector])/len(self.factor[detector]) * abs(float(np.linalg.det(cov)))**(1./4.)
      print emit, '\n'
      print "X Covariance:"
      print x_cov, "\n"
      x_emit = sum(self.factor[detector])/len(self.factor[detector]) * abs(float(np.linalg.det(x_cov)))**(1./2.)
      print x_emit, '\n'
      print "Y Covariance:"
      print y_cov, "\n"
      y_emit = sum(self.factor[detector])/len(self.factor[detector]) * abs(float(np.linalg.det(y_cov)))**(1./2.)
      print y_emit, '\n'

#########################################################################################
# Writes out the variable holding emittance terms to minimize what's in memory
#  def Clear_Data(self):
#    with open('emitt_data', 'a') as outfile:
#      for i in range(len(self.terms)):
#        json.dump(self.terms[i], outfile)

