import output_maus_analysis as _output
import numpy as np
import math

class Analysis(object):
  def __init__ (self):
    self.o_sel = _output.Output("Selection")
    _output.Message("INITALIZING EMITTANCE", exc=True)
    self.terms = []
    self.data = []
    self.factor = {"upstream":{1:[], 2:[], 3:[], 4:[], 5:[]}, \
                   "downstream":{1:[], 2:[], 3:[], 4:[], 5:[]}}

  def Process (self, data, cut_list):
    _output.Message("Collecting Emittance Data")
    if cut_list["Pass"]:
      event = {"upstream":{1:[], 2:[], 3:[], 4:[], 5:[]}, \
               "downstream":{1:[], 2:[], 3:[], 4:[], 5:[]}}
      for detector in data:
        tracks = data[detector]
        for i in range(len(tracks)):
          for j in range(len(tracks[i]["track_points"])):
            track_pt = tracks[i]["track_points"][j]
            if track_pt["plane"] == 0:
              station = track_pt["station"]
              _output.Message("Found a {} point".format(detector))
              event[detector][station] = [track_pt["x_pos"], track_pt["x_mom"], \
                                          track_pt["y_pos"], track_pt["y_mom"]]
              self.factor[detector][station].append(track_pt["z_mom"]/105.658)
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
    eval   = {"upstream": {1:[], 2:[], 3:[], 4:[], 5:[]}, \
              "downstream": {1:[], 2:[], 3:[], 4:[], 5:[]}}
    x_eval = {"upstream": {1:[], 2:[], 3:[], 4:[], 5:[]}, \
              "downstream": {1:[], 2:[], 3:[], 4:[], 5:[]}}
    y_eval = {"upstream": {1:[], 2:[], 3:[], 4:[], 5:[]}, \
              "downstream": {1:[], 2:[], 3:[], 4:[], 5:[]}}
    beam = {"upstream":[], "downstream":[]}; position = 0
    for event in self.terms:
      for detector in event:
        for station in event[detector]:
          eval[detector][station].append(event[detector][station])
          x_eval[detector][station].append([event[detector][station][0], event[detector][station][1]])
          y_eval[detector][station].append([event[detector][station][2], event[detector][station][3]])
    for detector in eval:
      for station in eval[detector]:

        # #######################################################################################################################
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
        #   first   += xpos_var[i] * xpos_var[i] / (size**2)
        #   second  += xpos_var[i] * xmom_var[i] / (size**2)
        #   third   += xpos_var[i] * ypos_var[i] / (size**2)
        #   fourth  += xpos_var[i] * ymom_var[i] / (size**2)
        #   fifth   += xmom_var[i] * xmom_var[i] / (size**2)
        #   sixth   += xmom_var[i] * ypos_var[i] / (size**2)
        #   seventh += xmom_var[i] * ymom_var[i] / (size**2)
        #   eighth  += ypos_var[i] * ypos_var[i] / (size**2)
        #   ninth   += ypos_var[i] * ymom_var[i] / (size**2)
        #   tenth   += ymom_var[i] * ymom_var[i] / (size**2)
        #
        # print (first)
        # print (second)
        # print (third)
        # print (fourth)
        # print (fifth)
        # print (sixth)
        # print (seventh)
        # print (eighth)
        # print (ninth)
        # print (tenth), "\n"
        #
        # #######################################################################################################################

        if detector == "upstream":
          if station == 1:
            position = 15060
          elif station == 2:
            position = 14860
          elif station == 3:
            position = 14610
          elif station == 4:
            position = 14310
          elif station == 5:
            position = 13960
        if detector == "downstream":
          if station == 1:
            position = 18850
          elif station == 2:
            position = 19050
          elif station == 3:
            position = 19300
          elif station == 4:
            position = 19600
          elif station == 5:
            position = 19950
        data   = np.array(eval[detector][station]).T
        x_data = np.array(x_eval[detector][station]).T
        y_data = np.array(y_eval[detector][station]).T
        cov    = np.cov(data)
        x_cov  = np.cov(x_data)
        y_cov  = np.cov(y_data)
        print detector.capitalize(), " Detector"
        print "Station ", station
        print "\nX Y Covariance:"
        print cov, "\n"
        emit   = sum(self.factor[detector])/len(self.factor[detector]) * abs(float(np.linalg.det(cov)))**(1./4.) / 105.658
        print emit, '\n'
        beam[detector].append([position, emit])
        self.o_sel.Graph("Emittance", "Emittance", position, emit)
        print "X Covariance:"
        print x_cov, "\n"
        x_emit = sum(self.factor[detector])/len(self.factor[detector]) * abs(float(np.linalg.det(x_cov)))**(1./2.) / 105.658
        print x_emit, '\n'
        print "Y Covariance:"
        print y_cov, "\n"
        y_emit = sum(self.factor[detector])/len(self.factor[detector]) * abs(float(np.linalg.det(y_cov)))**(1./2.) / 105.658
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

        # #################################################################################
        # x = []; px = []; y = []; py = []
        # for i in range(len(x_eval[detector])):
        #   x.append(x_eval[detector][i][0])
        #   px.append(x_eval[detector][i][1])
        #   y.append(y_eval[detector][i][0])
        #   py.append(y_eval[detector][i][1])
        # ave_pos_x = sum(x)  / len(x)
        # ave_mom_x = sum(px) / len(px)
        # ave_pos_y = sum(y)  / len(y)
        # ave_mom_y = sum(py) / len(py)
        # for i in range(len(x)):
        #   x[i]  = x[i]  - ave_pos_x
        #   px[i] = px[i] - ave_mom_x
        #   y[i]  = y[i]  - ave_pos_y
        #   py[i] = py[i] - ave_mom_y
        # x_fit = np.polyfit(x, px, 1, full=True)
        # y_fit = np.polyfit(y, py, 1, full=True)
        # theta_x = math.atan(x_fit[0][0])
        # theta_y = math.atan(y_fit[0][0])
        # for i in range(len(x_eval[detector])):
        #   x[i]  = math.cos(theta_x) * x[i]  + math.sin(theta_x)  * px[i]
        #   px[i] = math.cos(theta_x) * px[i] - math.sin(theta_x)  * x[i]
        #   y[i]  = math.cos(theta_y) * y[i]  + math.sin(theta_y)  * py[i]
        #   py[i] = math.cos(theta_y) * py[i] - math.sin(theta_y)  * y[i]
        #
        #   name = "Mod_Emit_X_{}".format(detector)
        #   title = "Modified Emittance X {}".format(detector)
        #   self.o_sel.Fill(name, title, x[i], px[i], \
        #                   200, -200, 200, 200, -200, 200)
        #   name = "Mod_Emit_Y_{}".format(detector)
        #   title = "Modified Emittance Y {}".format(detector)
        #   self.o_sel.Fill(name, title, y[i], py[i], \
        #                   200, -200, 200, 200, -200, 200)
        #
        #   name = "Mod_Pos_X_{}".format(detector)
        #   title = "Modified Position X {}".format(detector)
        #   self.o_sel.Fill(name, title, x[i], 200, -200, 200)
        #   name = "Mod_Pos_Y_{}".format(detector)
        #   title = "Modified Position Y {}".format(detector)
        #   self.o_sel.Fill(name, title, y[i], 200, -200, 200)
        #   name = "Mod_Mom_X_{}".format(detector)
        #   title = "Modified Momentum X {}".format(detector)
        #   self.o_sel.Fill(name, title, px[i], 200, -200, 200)
        #   name = "Mod_Mom_Y_{}".format(detector)
        #   title = "Modified Momentum Y {}".format(detector)
        #   self.o_sel.Fill(name, title, py[i], 200, -200, 200)
        #
        # print "Modified Fits:\n"
        # self.o_sel.Gaus_Fit(-200, 200, "Mod_Pos_X_{}".format(detector))
        # self.o_sel.Gaus_Fit(-200, 200, "Mod_Pos_Y_{}".format(detector))
        # self.o_sel.Gaus_Fit(-200, 200, "Mod_Mom_X_{}".format(detector))
        # self.o_sel.Gaus_Fit(-200, 200, "Mod_Mom_Y_{}".format(detector))
        #
        # y1 = x_fit[0][1] + x_fit[0][0] *  150.0
        # y2 = x_fit[0][1] + x_fit[0][0] * -150.0
        # print " 150, ", y1
        # print "-150, ", y2
        # #################################################################################