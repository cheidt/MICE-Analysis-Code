from config_maus_analysis import analysis_config as _config
import output_maus_analysis as _output
import numpy as np
import copy
import math as m
#########################################################################################
  # 
class ST_Alignment(object):
  def __init__ (self):
    _output.Message("INITALIZING ALIGNMENT", exc=True)
    self.o_res   = _output.Output("station_alignment")
    self.o_prof  = _output.Output("beam_profile")
    self.spaces  = {"upstream":{}, \
                    "downstream":{}}
    self.transit = {"upstream":{}, \
                    "downstream":{}}
    self.count   = 0
    self.position = {"upstream":{1:{"x":0, "y":0, "z":0, \
                                    "theta":0, "psi":0, "phi":0}, \
                                 2:{"x":0, "y":0, "z":0, \
                                    "theta":0, "psi":0, "phi":0}, \
                                 3:{"x":0, "y":0, "z":0, \
                                    "theta":0, "psi":0, "phi":0}, \
                                 4:{"x":0, "y":0, "z":0, \
                                    "theta":0, "psi":0, "phi":0}, \
                                 5:{"x":0, "y":0, "z":0, \
                                    "theta":0, "psi":0, "phi":0}}, \
                     "downstream":{1:{"x":0, "y":0, "z":0, \
                                      "theta":0, "psi":0, "phi":0}, \
                                   2:{"x":0, "y":0, "z":0, \
                                      "theta":0, "psi":0, "phi":0}, \
                                   3:{"x":0, "y":0, "z":0, \
                                      "theta":0, "psi":0, "phi":0}, \
                                   4:{"x":0, "y":0, "z":0, \
                                      "theta":0, "psi":0, "phi":0}, \
                                   5:{"x":0, "y":0, "z":0, \
                                      "theta":0, "psi":0, "phi":0}}}

#########################################################################################
  # 
  def StS_Collect_Space_Points(self, tracks):
    for detector in tracks:
      size = len(tracks[detector])
      if not size == 1:
        if size == 0:
          continue
        else:
          _output.Message("Number of tracks in the ", detector, " ", size)
          pop_list = []
          for track in range(len(tracks[detector])):
            _output.Message("Track ", track, " has ", tracks[detector]\
                            [track]["triplets"], " data points")
            if int(tracks[detector][track]["triplets"]) <= \
               int(_config["min_trip"]):
              pop_list.append(track)
          pop_size = len(pop_list)
          if not pop_size == 0:
            pop_list.sort(reverse=True)
            for i in range(pop_size):
              tracks[detector].pop(pop_list[i])

      size = len(tracks[detector])
      if not size == 1:
        if size == 0:
          _output.Message("No acceptable tracks")
          continue
        if size > 1:
          _output.Message("Too many tracks")
          continue

      if tracks[detector][0]["triplets"] < _config["req_trip"]:
        continue

      for station in tracks[detector][0]["seeds"][0][detector]:
        if not len(tracks[detector][0]["seeds"][0][detector][station]) == 1:
          _output.Message("Phantom space point in ",detector," ",station, \
                                                    exc=True)

      for seeds in tracks[detector][0]["seeds"]:
        temp = np.array([(0,0,0,0,0),(0,0,0,0,0),(0,0,0,0,0), \
                         (0,0,0,0,0),(0,0,0,0,0)],
                  dtype=[('x_pos','f4'),('y_pos','f4'),('z_pos','f4'), \
                         ('x_exp','f4'),('y_exp','f4')])
        use_event = True
        for station in seeds[detector]:
          try:
            x_pos   = seeds[detector][station][0]["x_glob_pos"]
            y_pos   = seeds[detector][station][0]["y_glob_pos"]
            z_pos   = seeds[detector][station][0]["z_glob_pos"]
            i       = station - 1
            temp[i] = (x_pos, y_pos, z_pos, 0, 0)
          except IndexError:
            use_event = False
            continue

        if use_event:
          if len(self.spaces[detector]) > 0:
            self.spaces[detector] = np.vstack((self.spaces[detector],temp))
          else:
            self.spaces[detector] = temp

#########################################################################################
  # 
  def Station_Alignment(self):
    self.transit = copy.deepcopy(self.spaces)
    temp         = copy.deepcopy(self.spaces)
    loop = 0

    while loop < 1:
      loop += 1
      _output.Message("Loop over data number ", loop, exc=True)
      for detector in self.spaces:
        delete_list = self.Fiducial_Cut(temp[detector])
        temp[detector] = np.delete(temp[detector], delete_list, 0)
        for cut_station in range(5):
          delete_list = self.Draw_Line(temp[detector], cut_station)
          temp[detector] = np.delete(temp[detector], delete_list, 0)
        for cut_station in range (5):
          print "\n", detector, " Tracker"
          print "Station ", cut_station
          self.Find_Coefficents(temp[detector],cut_station)
      #print self.transit
          #x_res        = x_exp[event] - space["x_pos"]
          #y_res        = y_exp[event] - space["y_pos"]
          #res[event]   = (y_res**2 + x_res**2)**0.5

          #name  = "st_st_residual"
          #title = "Four Station Residual"
          #self.o_res.Fill(name, title, x_res, y_res, \
          #                500, -50 , 50, 500, -50 , 50, \
          #                detector=detector, station=cut_station, \
          #                iteration=loop)
          #name  = "st_st_profile"
          #title = "Beam Profile in Space Points"
          #self.o_prof.Fill(name, title, space["x_pos"], space["y_pos"], \
          #                 500, -250 , 250, 500, -250 , 250, \
          #                 detector=detector, station=cut_station, \
          #                 iteration=loop)

        #chi_sq[loop][detector][cut_station] = self.Chi_Squared(res)
        #fit[detector][cut_station] = self.Find_Coefficents(\
        #                             self.transit[detector], x_exp, y_exp, \
        #                                  cut_station)

      #self.Move_Points(fit)

    #self.o_res.Write()
    #self.o_prof.Write()

#########################################################################################
  #  Cuts on a certain area within the trackers, requires all parts of the particle
  #    path to pass within that area.
  def Fiducial_Cut(self, seeds):
    delete_list = []
    distance = (seeds['x_pos']**2 + seeds['y_pos']**2)**0.5
    for i in range(len(distance)):
      for j in range(len(distance[i])):
        if distance[i][j] < 100.0:
          delete_list.append(i)
          break
    return delete_list


#########################################################################################
  #  Draws a line between any four space points and determines where in the fifth plane
  #    we would expect to find the space point
  def Draw_Line(self, seeds, cut):
    cut_lst = [0, 1, 2, 3, 4]
    delete_list = []
    x_array = np.delete(seeds['x_pos'],cut,1)
    y_array = np.delete(seeds['y_pos'],cut,1)
    z_array = np.delete(seeds['z_pos'],cut,1)
    cut_lst.pop(cut)
    z_value = np.delete(seeds['z_pos'],cut_lst,1)

    for i in range(len(x_array)):
      Z = np.vstack([z_array[i], np.ones(len(z_array[i]))]).T
      x_result = np.linalg.lstsq(Z, x_array[i])
      seeds[i][cut]["x_exp"] = z_value[i]*x_result[0][0] + x_result[0][1]
      y_result = np.linalg.lstsq(Z, y_array[i])
      seeds[i][cut]["y_exp"] = z_value[i]*y_result[0][0] + y_result[0][1]

      if x_result[1] > 3 or y_result[1] > 3:
        delete_list.append(i)
        
    return delete_list

#########################################################################################
  #  Takes expected and actual values and returns the rotation matrix that transforms
  #    one to the other.
  def Find_Coefficents(self, seeds, cut):
    station = [0, 1, 2, 3, 4]
    station.pop(cut)
    x_value  = np.delete(seeds['x_pos'],station,1)
    y_value  = np.delete(seeds['y_pos'],station,1)
    z_value  = np.delete(seeds['z_pos'],station,1)
    x_expect = np.delete(seeds['x_exp'],station,1) # - x_value
    y_expect = np.delete(seeds['y_exp'],station,1) # - y_value

    X = np.hstack([x_value, y_value, np.ones((len(z_value),1))])
    Y = np.hstack([x_value, y_value, np.ones((len(z_value),1))])
    Z = np.hstack([y_value, np.ones((len(z_value),1))])

    print "Number of Points - ", len(x_value),"\n"

    print "Number of dimensions in set"
    print np.shape(X)
    print np.shape(x_expect),"\n"

    print "Average Residual"
    print "X - ", np.average(x_expect - x_value)
    print "Y - ", np.average(y_expect - y_value), "\n"

    use_individual = 0

    if (use_individual == 1):
      print "Using Indiv"
      for t in range(len(X)):
        print "Event -",t

        oneD_X  = np.expand_dims(X[t], 0)
        oneD_Y  = np.expand_dims(Y[t], 0)
        oneD_xe = np.expand_dims(x_expect[t], 0)
        oneD_ye = np.expand_dims(y_expect[t], 0)

        print "Data points"
        print "X Values"
        print oneD_X
        print oneD_xe, "\n"
        print "Y Values"
        print oneD_Y
        print oneD_ye, "\n"

        a12, xt = np.linalg.lstsq(oneD_X, oneD_xe)[0]
        a21, yt = np.linalg.lstsq(oneD_Y, oneD_ye)[0]

        coeff = np.array([[1, a12, 0, xt], [a21, 1, 0, yt]])
        print "Rotation Matrix\n",coeff, "\n"

        test_x = 1 * x_value[t] + a12 * y_value[t] + 0 * z_value[t] + xt
        print "Expected X value:   ", x_value[t]
        print "Calculated X value: ", test_x

        test_y = a21 * x_value[t] + 1 * y_value[t] + 0 * z_value[t] + yt
        print "Expected Y value:   ", y_value[t]
        print "Calculated Y value: ", test_y, "\n"
      
    if (use_individual == 0):
      a11, a12, xt = np.linalg.lstsq(X, x_expect)[0]
      a21, a22, yt = np.linalg.lstsq(Y, y_expect)[0]
      #a31, a32, zt = np.linalg.lstsq(Z, x_value,1)))[0]

      calc_a21 = -a23 * a13

      coeff = np.array([[a11[0], a12, 0, xt], [a21, a22, 0, yt], [a31, a32, 0, zt]])
      print "Rotation Matrix\n",coeff, "\n"
      
#      for t in range(len(x_value)):
#        test_x = 1 * x_value[t] + a12 * y_value[t] + a13 * z_value[t] + xt
#        print "Expected X value:   ", x_value[t]
#        print "Calculated X value: ", test_x

#        test_y = a21 * x_value[t] + 1 * y_value[t] + a23 * z_value[t] + yt
#        print "Expected Y value:   ", y_value[t]
#        print "Calculated Y value: ", test_y

#########################################################################################
  # 
  def Chi_Squared(self, x):
    mean   = sum(x.values())/len(x)
    n      = len(x)
    chi_sq = 0
    for i in x:
      chi_sq += ((x[i]-mean)**2)/(mean*n)
    return chi_sq

#########################################################################################
  # 
  def Move_Points(self, fit):
    step = .25
    
    for detector in self.transit:
      for station in range (1,6):
        x     =  fit[detector][station].item(3)
        y     =  fit[detector][station].item(7)
        z     =  0
        phi   =  fit[detector][station].item(1)
        psi   =  fit[detector][station].item(6)
        theta = -fit[detector][station].item(2)

        sx  = self.position[detector][station]["x"]     = \
              step * (x - self.position[detector][station]["x"]) + \
              self.position[detector][station]["x"]
        sy  = self.position[detector][station]["y"]     = \
              step * (y - self.position[detector][station]["y"]) + \
              self.position[detector][station]["y"]
        sz  = self.position[detector][station]["z"]     = \
              step * (z - self.position[detector][station]["z"]) + \
              self.position[detector][station]["z"]
        sph = self.position[detector][station]["phi"]   = \
              (phi - self.position[detector][station]["phi"]) + \
              self.position[detector][station]["phi"]
        sps = self.position[detector][station]["psi"]   = \
              (psi - self.position[detector][station]["psi"]) + \
              self.position[detector][station]["psi"]
        sth = self.position[detector][station]["theta"] = \
              (theta - self.position[detector][station]["theta"]) + \
              self.position[detector][station]["theta"]
     
        step_fit = np.matrix([[m.cos(sth)*m.cos(sph), m.cos(sth)*m.sin(sph), \
                              -m.sin(sth), sx], \
                              [m.sin(sps)*m.sin(sth)*m.cos(sph) - m.cos(sps)*m.sin(sph), \
                               m.sin(sps)*m.sin(sth)*m.sin(sph) + m.cos(sps)*m.cos(sph), \
                               m.cos(sth)*m.sin(sps), sy], \
                              [m.cos(sps)*m.sin(sth)*m.cos(sph) + m.sin(sps)*m.sin(sph), \
                               m.cos(sps)*m.sin(sth)*m.sin(sph) - m.sin(sps)*m.cos(sph), \
                               m.cos(sth)*m.cos(sps), sz], \
                              [0, 0, 0, 1]])

        for event in self.transit[detector]:
          point = self.transit[detector][event][station]
          tran = np.matrix([[point["x_pos"]],[point["y_pos"]], \
                            [point["z_pos"]],[1]])
          temp = step_fit*tran

          point["x_pos"] += step * (temp.item(0) - point["x_pos"])
          point["y_pos"] += step * (temp.item(1) - point["y_pos"])
          point["z_pos"] += step * (temp.item(2) - point["z_pos"])
