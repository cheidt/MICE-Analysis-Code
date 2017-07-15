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

  def Process (self, data, cut_list):
    _output.Message("Collecting Emittance Data")
    if cut_list["Pass"]:
      event = {"upstream":[], "downstream":[]}
      for detector in data:
        tracks = data[detector]
        for i in range(len(tracks)):
          for j in range(len(tracks[i]["track_points"])):
            track_pt = tracks[i]["track_points"][j]
            if track_pt["station"] == 5 and track_pt["plane"] == 0 and \
               track_pt["x_mom"]**2 > .000001  and track_pt["y_mom"]**2 > .000001:
              _output.Message("Found a {} point".format(detector))
              event[detector] = [track_pt["x_pos"], track_pt["x_mom"], \
                                 track_pt["y_pos"], track_pt["y_mom"]]
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
      if self.terms[-1]["upstream"] == []:
        print "upstream"
        print "event: ", len(self.terms)
        print self.terms[-1]["upstream"], "\n\n"
      if self.terms[-1]["downstream"] == []:
        print "downstream"
        print "event: ", len(self.terms)
        print self.terms[-1]["downstream"], "\n\n"

  def Write(self):
    self.o_sel.Write()

  #########################################################################################
  # Takes in emittance terms and outputs system emittance.  Separates out cut booleans
  # and uses those to determine what makes it to final calculation.  Solves for
  # 4D and the two 2D emittances.
  #########################################################################################
  def Emittance(self):
    eval   = {"upstream": [], "downstream": []}
    x_eval = {"upstream": [], "downstream": []}
    y_eval = {"upstream": [], "downstream": []}
    for event in self.terms:
      for detector in event:
        eval[detector].append(event[detector])
        # print event[detector]
        x_eval[detector].append([event[detector][0], event[detector][1]])
        y_eval[detector].append([event[detector][2], event[detector][3]])
    for detector in eval:
########################################################################################################################
      # print "Check"
      # xpos = []; xmom = []; xpos_sqr = []; xmom_sqr = []
      # test1 = []; test2 = []; test3 = []; test4 = []; test5 = []
      # ypos = []; ymom = []; ypos_sqr = []; ymom_sqr = []
      # size = len(x_eval[detector]) - 1
      # for i in range(size):
      #   xpos.append(x_eval[detector][i][0])
      #   xpos_sqr.append(x_eval[detector][i][0]**2)
      #   xmom.append(x_eval[detector][i][1])
      #   xmom_sqr.append(x_eval[detector][i][1]**2)
      #   ypos.append(y_eval[detector][i][0])
      #   ypos_sqr.append(y_eval[detector][i][0]**2)
      #   ymom.append(y_eval[detector][i][1])
      #   ymom_sqr.append(y_eval[detector][i][1]**2)
      #   test2.append(x_eval[detector][i][0]*y_eval[detector][i][0])
      #   test3.append(x_eval[detector][i][0]*y_eval[detector][i][1])
      #   test1.append(x_eval[detector][i][0]*x_eval[detector][i][1])
      #   test4.append(x_eval[detector][i][1]*y_eval[detector][i][0])
      #   test5.append(x_eval[detector][i][1]*y_eval[detector][i][1])
      #
      # xave_pos = sum(xpos) / size
      # xave_mom = sum(xmom) / size
      # xave_pos_sqr = sum(xpos_sqr) / size
      # xave_mom_sqr = sum(xmom_sqr) / size
      # yave_pos = sum(ypos) / size
      # yave_mom = sum(ymom) / size
      # yave_pos_sqr = sum(ypos_sqr) / size
      # yave_mom_sqr = sum(ymom_sqr) / size
      # ave_test1 = sum(test1) / size
      # ave_test2 = sum(test2) / size
      # ave_test3 = sum(test3) / size
      # ave_test4 = sum(test4) / size
      # ave_test5 = sum(test5) / size
      #
      # print "X Pos: ", xave_pos, " / ", xave_pos_sqr
      # print "X Mom: ", xave_mom, " / ", xave_mom_sqr
      # print "Y Pos: ", yave_pos, " / ", yave_pos_sqr
      # print "Y Mom: ", yave_mom, " / ", yave_mom_sqr
      #
      # print "Test1: ", ave_test1
      # print "Test2: ", ave_test2
      # print "Test3: ", ave_test3
      # print "Test4: ", ave_test4
      # print "Test5: ", ave_test5, "\n"
      #
      # first = abs((xave_pos**2 - xave_pos_sqr))**(.5)
      # fourth = abs((xave_mom**2 - xave_mom_sqr))**(.5)
      # ninth = abs((yave_pos**2 - yave_pos_sqr))**(.5)
      # sixteenth = abs((yave_mom**2 - yave_mom_sqr))**(.5)
      #
      # print "First: ", first
      # print "Fourth: ", fourth
      # print "Ninth: ", ninth
      # print "Sixteenth: ", sixteenth, "\n"
      #
      # xpos_var = [i - xave_pos for i in xpos]
      # xmom_var = [i - xave_mom for i in xmom]
      # ypos_var = [i - yave_pos for i in ypos]
      # ymom_var = [i - yave_mom for i in ymom]
      #
      # first = 0; second = 0; third = 0; fourth = 0
      # fifth = 0; sixth = 0; seventh = 0; eighth = 0
      # ninth = 0; tenth = 0
      # for i in range(len(xpos_var)):
      #   first   += xpos_var[i] * xpos_var[i]
      #   second  += xpos_var[i] * xmom_var[i]
      #   third   += xpos_var[i] * ypos_var[i]
      #   fourth  += xpos_var[i] * ymom_var[i]
      #   fifth   += xmom_var[i] * xmom_var[i]
      #   sixth   += xmom_var[i] * ypos_var[i]
      #   seventh += xmom_var[i] * ymom_var[i]
      #   eighth  += ypos_var[i] * ypos_var[i]
      #   ninth   += ypos_var[i] * ymom_var[i]
      #   tenth   += ymom_var[i] * ymom_var[i]
      #
      # print (first / (size))
      # print (second / (size))
      # print (third / (size))
      # print (fourth / (size))
      # print (fifth / (size))
      # print (sixth / (size))
      # print (seventh / (size))
      # print (eighth / (size))
      # print (ninth / (size))
      # print (tenth / (size)), "\n"

########################################################################################################################

      data   = np.array(eval[detector]).T
      x_data = np.array(x_eval[detector]).T
      y_data = np.array(y_eval[detector]).T
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

########################################################################################################################
      # print "Experiment:"
      # print "\nX Y Covariance:"
      # print cov, "\n"
      # emit = sum(self.factor[detector]) / len(self.factor[detector]) * abs(float(np.linalg.det(cov))) ** (1. / 8.)
      # print emit, '\n'
      # print "X Covariance:"
      # print x_cov, "\n"
      # x_emit = sum(self.factor[detector]) / len(self.factor[detector]) * abs(float(np.linalg.det(x_cov))) ** (1. / 4.)
      # print x_emit, '\n'
      # print "Y Covariance:"
      # print y_cov, "\n"
      # y_emit = sum(self.factor[detector]) / len(self.factor[detector]) * abs(float(np.linalg.det(y_cov))) ** (1. / 4.)
      # print y_emit, '\n'
      #
      # print "Factor: ", sum(self.factor[detector]) / len(self.factor[detector])
########################################################################################################################

      xpos = []; ypos = []; xmom = []; ymom = []
      for detector in x_eval:
        for i in x_eval[detector]:
          xpos.append(x_eval[detector][i][0])
          x_fit = np.polyfit(x_eval[detector][1], x_eval[detector][0], 1, full=True)