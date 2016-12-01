from config_maus_analysis import analysis_config as _config
import output_maus_analysis as _output
import numpy as np
import copy
import math as m


#########################################################################################
#
class ST_Alignment(object):
  def __init__(self):
    _output.Message("INITIALIZING ALIGNMENT", exc = True)
    self.o_res = _output.Output("station_alignment")
    self.o_prof = _output.Output("beam_profile")
    self.spaces = {"upstream"  : {}, \
                   "downstream": {}}

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
            _output.Message("Track ", track, " has ", tracks[detector] \
              [track]["triplets"], " data points")
            if int(tracks[detector][track]["triplets"]) <= \
                    int(_config["min_trip"]):
              pop_list.append(track)
          pop_size = len(pop_list)
          if not pop_size == 0:
            pop_list.sort(reverse = True)
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
          _output.Message(detector, "Phantom space point in ", " ", station, \
                          exc = True)

      for seeds in tracks[detector][0]["seeds"]:
        temp = np.array([(0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0), \
                         (0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0)],
                        dtype = [('x_pos', 'f4'), ('y_pos', 'f4'), ('z_pos', 'f4'), \
                                 ('x_exp', 'f4'), ('y_exp', 'f4'), ('z_exp', 'f4')])
        use_event = True
        for station in seeds[detector]:
          try:
            x_pos = seeds[detector][station][0]["x_pos"]
            y_pos = seeds[detector][station][0]["y_pos"]
            z_pos = seeds[detector][station][0]["z_pos"]
            i = station - 1
            temp[i] = (x_pos, y_pos, z_pos, 0, 0, 0)
          except IndexError:
            use_event = False
            continue

        if use_event:
          if len(self.spaces[detector]) > 0:
            self.spaces[detector] = np.vstack((self.spaces[detector], temp))
          else:
            self.spaces[detector] = temp
    pass

  #########################################################################################
  # 
  def Station_Alignment(self):
    self.transit = copy.deepcopy(self.spaces)
    transform = {"upstream": {0: {}, 1: {}, 2: {}, 3: {}, 4: {}}, \
               "downstream": {0: {}, 1: {}, 2: {}, 3: {}, 4: {}}}
    loop = 0

    while loop < 2:
      temp = copy.deepcopy(self.transit)
      loop += 1
      _output.Message("Loop over data number ", loop, exc = True)
      for detector in self.spaces:
        delete_list = self.Fiducial_Cut(temp[detector])
        temp[detector] = np.delete(temp[detector], delete_list, 0)
        for cut_station in range(5):
          delete_list = self.Find_Expected_Values(temp[detector], cut_station)
          temp[detector] = np.delete(temp[detector], delete_list, 0)
        for cut_station in range(5):
          result = self.Find_Coefficents(temp[detector], cut_station)
          transform[detector][cut_station] = result
      print "Final Transform:"
      print transform
      x = {"g":5, "j": 6}
      transform
      self.Move_Points(transform)
      # print self.transit
      # x_res        = x_exp[event] - space["x_pos"]
      # y_res        = y_exp[event] - space["y_pos"]
      # res[event]   = (y_res**2 + x_res**2)**0.5

      # name  = "st_st_residual"
      # title = "Four Station Residual"
      # self.o_res.Fill(name, title, x_res, y_res, \
      #                500, -50 , 50, 500, -50 , 50, \
      #                detector=detector, station=cut_station, \
      #                iteration=loop)
      # name  = "st_st_profile"
      # title = "Beam Profile in Space Points"
      # self.o_prof.Fill(name, title, space["x_pos"], space["y_pos"], \
      #                 500, -250 , 250, 500, -250 , 250, \
      #                 detector=detector, station=cut_station, \
      #                 iteration=loop)

      # chi_sq[loop][detector][cut_station] = self.Chi_Squared(res)
      # fit[detector][cut_station] = self.Find_Coefficents(\
      #                             self.transit[detector], x_exp, y_exp, \
      #                                  cut_station)

    # self.o_res.Write()
    # self.o_prof.Write()
    pass

  #########################################################################################
  #  Cuts on a certain area within the trackers, requires all parts of the particle
  #    path to pass within that area.
  def Fiducial_Cut( self, seeds ):
    delete_list = []
    distance = (seeds['x_pos'] ** 2 + seeds['y_pos'] ** 2) ** 0.5
    for i in range(len(distance)):
      for j in range(len(distance[i])):
        if distance[i][j] < 0.0:
          delete_list.append(i)
          break
    return delete_list

  #########################################################################################
  #  Draws a line between any four space points and determines where in the fifth plane
  #    we would expect to find the space point
  def Find_Expected_Values(self, seeds, cut):
    cut_list = [0, 1, 2, 3, 4]
    cut_list.pop(cut)
    delete_list = []
    x_array = np.delete(seeds['x_pos'], cut, 1)
    y_array = np.delete(seeds['y_pos'], cut, 1)
    z_array = np.delete(seeds['z_pos'], cut, 1)
    x_value = np.delete(seeds['x_pos'], cut_list, 1)
    y_value = np.delete(seeds['y_pos'], cut_list, 1)
    z_value = np.delete(seeds['z_pos'], cut_list, 1)

    for i in range(len(x_array)):
      ZX = np.vstack([z_array[i], x_array[i], np.ones(len(z_array[i]))]).T
      ZY = np.vstack([z_array[i], y_array[i], np.ones(len(z_array[i]))]).T
      XY = np.vstack([x_array[i], y_array[i], np.ones(len(z_array[i]))]).T
      x_result = np.linalg.lstsq(ZY, x_array[i])
      y_result = np.linalg.lstsq(ZX, y_array[i])
      z_result = np.linalg.lstsq(XY, z_array[i])

      seeds[i][cut]["x_exp"] = z_value[i] * x_result[0][0] + y_value[i] * x_result[0][1] + \
                               x_result[0][2]
      seeds[i][cut]["y_exp"] = z_value[i] * y_result[0][0] + x_value[i] * y_result[0][1] + \
                               y_result[0][2]
      seeds[i][cut]["z_exp"] = x_value[i] * z_result[0][0] + y_value[i] * z_result[0][1] + \
                               z_result[0][2]

      if x_result[1] > 10 or y_result[1] > 10:
        delete_list.append(i)

    return delete_list

  #########################################################################################
  #  Takes expected and actual values and returns the rotation matrix that transforms
  #    one to the other.
  def Find_Coefficents( self, seeds, cut ):
    station = [0, 1, 2, 3, 4]
    station.pop(cut)
    x_value = np.delete(seeds['x_pos'], station, 1)
    y_value = np.delete(seeds['y_pos'], station, 1)
    z_value = np.delete(seeds['z_pos'], station, 1)
    x_expect = np.delete(seeds['x_exp'], station, 1)  # - x_value
    y_expect = np.delete(seeds['y_exp'], station, 1)  # - y_value
    z_expect = np.delete(seeds['z_exp'], station, 1)

    X = np.hstack([x_value, y_value, z_value, np.ones((len(z_value), 1))])
    Y = np.hstack([x_value, y_value, z_value, np.ones((len(z_value), 1))])
    Z = np.hstack([x_value, y_value, z_value, np.ones((len(z_value), 1))])

    print "Number of Points - ", len(x_value), "\n"

    print "Average Residual"
    print "X - ", np.average(x_expect - x_value)
    print "Y - ", np.average(y_expect - y_value)
    print "Z - ", np.average(z_expect - z_value), "\n"

    a11, a12, a13, xt = np.linalg.lstsq(X, x_expect)[0]
    a21, a22, a23, yt = np.linalg.lstsq(Y, y_expect)[0]
    a31, a32, a33, zt = np.linalg.lstsq(Z, z_expect)[0]

    coeff = np.array([[a11[0], a12, a13, xt], [a21, a22, a23, yt], [a31, a32, a33, zt]])
    print "Rotation Matrix\n", coeff, "\n"

    theta = -a13[0]
    phi = a12[0]
    psi = a23[0]

    #print 'checking: First order expect'
    #print 'theta: ', theta
    #print 'phi:   ', phi
    #print 'psi:   ', psi
    #print 'checking: Second order expect'
    #print 'theta: ', a31 - (phi * psi)
    #print 'phi:   ', -a21 + (psi * theta)
    #print 'psi:   ', -a32 + (phi * theta)

    out = {"phi": phi, "psi": psi, "theta": theta, "x": xt[0], "y": yt[0], "z": zt[0]}

    #    for t in range(len(x_value)):
    #      test_x = a21 * x_value[t] + a12 * y_value[t] + a13 * z_value[t] + xt
    #      print "Expected X value:   ", x_value[t]
    #      print "Calculated X value: ", test_x
    #
    #      test_y = a21 * x_value[t] + a22 * y_value[t] + a23 * z_value[t] + yt
    #      print "Expected Y value:   ", y_value[t]
    #      print "Calculated Y value: ", test_y
    #
    #      test_z = a31 * x_value[t] + a32 * y_value[t] + a33 * z_value[t] + zt
    #      print "Expected Z value:   ", z_value[t]
    #      print "Calculated Z value: ", test_z

    return out

  #########################################################################################
  # 
  def Move_Points(self, transform):
    step = .10

    for detector in self.transit:
      for station in range(0, 5):
        sx  = transform[detector][station]['x'] * step
        sy  = transform[detector][station]['y'] * step
        sz  = transform[detector][station]['z'] * step
        sph = transform[detector][station]['phi'] * step
        sps = transform[detector][station]['psi'] * step
        sth = transform[detector][station]['theta'] * step

        step_fit = np.matrix([[m.cos(sth) * m.cos(sph), m.cos(sth) * m.sin(sph), \
                               -m.sin(sth), sx], \
                              [m.sin(sps) * m.sin(sth) * m.cos(sph) - m.cos(sps) * m.sin(sph), \
                               m.sin(sps) * m.sin(sth) * m.sin(sph) + m.cos(sps) * m.cos(sph), \
                               m.cos(sth) * m.sin(sps), sy], \
                              [m.cos(sps) * m.sin(sth) * m.cos(sph) + m.sin(sps) * m.sin(sph), \
                               m.cos(sps) * m.sin(sth) * m.sin(sph) - m.sin(sps) * m.cos(sph), \
                               m.cos(sth) * m.cos(sps), sz], \
                              [0, 0, 0, 1]])

        print "Tracker ", detector
        print "Station ", station
        print step_fit

        for event in self.transit[detector]:
          point = event[station]
          tran = np.matrix([[point["x_pos"]], [point["y_pos"]], \
                            [point["z_pos"]], [1]])
          temp = step_fit * tran

          point["x_pos"] = temp.item(0)
          point["y_pos"] = temp.item(1)
          point["z_pos"] = temp.item(2)