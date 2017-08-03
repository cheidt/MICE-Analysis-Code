import output_maus_analysis as _output
import numpy as np
import math

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
            if track_pt["station"] == 5 and track_pt["plane"] == 0:
              _output.Message("Found a {} point".format(detector))
              event[detector] = [track_pt["x_pos"], track_pt["x_mom"], \
                                 track_pt["y_pos"], track_pt["y_mom"]]
              self.factor[detector].append(track_pt["z_mom"]/105.66)
              name = "Emit_X_{}".format(detector)
              title = "Emittance X {}".format(detector)
              self.o_sel.Fill(name, title, track_pt["x_pos"], track_pt["x_mom"], \
                              200, -200, 200, 200, -200, 200)
              name = "Emit_Y_{}".format(detector)
              title = "Emittance Y {}".format(detector)
              self.o_sel.Fill(name, title, track_pt["y_pos"], track_pt["y_mom"], \
                              200, -200, 200, 200, -200, 200)
          if len(tracks) > 1:
            print detector
            print "Why so many tracks?"
            raw_input("Press Enter to Exit")
          if event[detector] == []:
            print detector
            print event[detector]
            print cut_list
            print tracks[i]["track_points"]
            raw_input("Press Enter to Exit")
      self.terms.append(event)
      if self.terms[-1] == []:
        print "Problem writting event: ", len (self.terms)
        print self.terms[-1]


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
    i = 0
    for event in self.terms:
      for detector in event:
        eval[detector].append(event[detector])
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

      x = []; px = []; y = []; py = []
      for i in range(len(x_eval[detector])):
        x.append(x_eval[detector][i][0])
        px.append(x_eval[detector][i][1])
        y.append(y_eval[detector][i][0])
        py.append(y_eval[detector][i][1])
      ave_pos_x = sum(x)  / len(x)
      ave_mom_x = sum(px) / len(px)
      ave_pos_y = sum(y)  / len(y)
      ave_mom_y = sum(py) / len(py)
      for i in range(len(x)):
        x[i]  = x[i]  - ave_pos_x
        px[i] = px[i] - ave_mom_x
        y[i]  = y[i]  - ave_pos_y
        py[i] = py[i] - ave_mom_y
      x_fit = np.polyfit(x, px, 1, full=True)
      y_fit = np.polyfit(y, py, 1, full=True)
      theta_x = math.atan(x_fit[0][0])
      theta_y = math.atan(y_fit[0][0])
      for i in range(len(x_eval[detector])):
        x[i]  = math.cos(theta_x) * x[i]  + math.sin(theta_x)  * px[i]
        px[i] = math.cos(theta_x) * px[i] - math.sin(theta_x)  * x[i]
        y[i]  = math.cos(theta_y) * y[i]  + math.sin(theta_y)  * py[i]
        py[i] = math.cos(theta_y) * py[i] - math.sin(theta_y)  * y[i]

        name = "Mod_Emit_X_{}".format(detector)
        title = "Modified Emittance X {}".format(detector)
        self.o_sel.Fill(name, title, x[i], px[i], \
                        200, -200, 200, 200, -200, 200)
        name = "Mod_Emit_Y_{}".format(detector)
        title = "Modified Emittance Y {}".format(detector)
        self.o_sel.Fill(name, title, y[i], py[i], \
                        200, -200, 200, 200, -200, 200)

        name = "Mod_Pos_X_{}".format(detector)
        title = "Modified Position X {}".format(detector)
        self.o_sel.Fill(name, title, x[i], 200, -200, 200)
        name = "Mod_Pos_Y_{}".format(detector)
        title = "Modified Position Y {}".format(detector)
        self.o_sel.Fill(name, title, y[i], 200, -200, 200)
        name = "Mod_Mom_X_{}".format(detector)
        title = "Modified Momentum X {}".format(detector)
        self.o_sel.Fill(name, title, px[i], 200, -200, 200)
        name = "Mod_Mom_Y_{}".format(detector)
        title = "Modified Momentum Y {}".format(detector)
        self.o_sel.Fill(name, title, py[i], 200, -200, 200)

      print "Modified Fits:\n"
      self.o_sel.Gaus_Fit(-200, 200, "Mod_Pos_X_{}".format(detector))
      self.o_sel.Gaus_Fit(-200, 200, "Mod_Pos_Y_{}".format(detector))
      self.o_sel.Gaus_Fit(-200, 200, "Mod_Mom_X_{}".format(detector))
      self.o_sel.Gaus_Fit(-200, 200, "Mod_Mom_Y_{}".format(detector))

      y1 = x_fit[0][1] + x_fit[0][0] *  150.0
      y2 = x_fit[0][1] + x_fit[0][0] * -150.0
      print " 150, ", y1
      print "-150, ", y2