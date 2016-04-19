from config_maus_analysis import analysis_config as _config
import output_maus_analysis as _output
import numpy as np
import copy
import math as m
#########################################################################################
  # 
class ST_Alignment(object):
  def __init__ (self):
    print "INITALIZING ALIGNMENT"
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
        temp = np.array([(0,0,0),(0,0,0),(0,0,0),(0,0,0),(0,0,0)],
                  dtype=[('x_pos','f4'),('y_pos','f4'),('z_pos','f4'),])
        use_event = True
        for station in seeds[detector]:
          try:
            x_pos   = seeds[detector][station][0]["x_glob_pos"]
            y_pos   = seeds[detector][station][0]["y_glob_pos"]
            z_pos   = seeds[detector][station][0]["z_glob_pos"]
            i       = station - 1
            temp[i] = (x_pos, y_pos, z_pos)
          except IndexError:
            use_event = False
            print "Don't use event"
            continue

        if use_event:
          if len(self.spaces[detector]) > 0:
            self.spaces[detector] = np.vstack((self.spaces[detector],temp))
          else:
            self.spaces[detector] = temp
        else:
          print "Not Using event"

#########################################################################################
  # 
  def Station_Alignment(self):
    self.transit = copy.deepcopy(self.spaces)
    chi_sq = {}
    loop = 0
    while loop < 15:
      _output.Message("Loop over data number ", loop, exc=True)
      loop += 1
      chi_sq[loop] = {"upstream":{},
                      "downstream":{}}
      fit = {}
      for detector in self.spaces:
        fit[detector] = {}
        for cut_station in range(5):
          x_exp, y_exp, res = {}, {}, {}

#          space = self.transit[detector][event][cut_station]
#          temp = copy.deepcopy(self.transit[detector][event])
#          temp.pop(cut_station, None)
          line = self.Draw_Line(self.transit, cut_station)
          x_exp[event] = float(line["mx"]*space["z_pos"]+line["bx"])
          x_res        = x_exp[event] - space["x_pos"]
          y_exp[event] = float(line["my"]*space["z_pos"]+line["by"])
          y_res        = y_exp[event] - space["y_pos"]
          res[event]   = (y_res**2 + x_res**2)**0.5

          name  = "st_st_residual"
          title = "Four Station Residual"
          self.o_res.Fill(name, title, x_res, y_res, \
                          500, -50 , 50, 500, -50 , 50, \
                          detector=detector, station=cut_station, \
                          iteration=loop)
          name  = "st_st_profile"
          title = "Beam Profile in Space Points"
          self.o_prof.Fill(name, title, space["x_pos"], space["y_pos"], \
                           500, -250 , 250, 500, -250 , 250, \
                           detector=detector, station=cut_station, \
                           iteration=loop)

        #chi_sq[loop][detector][cut_station] = self.Chi_Squared(res)
        fit[detector][cut_station] = self.Find_Coefficents(\
                                     self.transit[detector], x_exp, y_exp, \
                                          cut_station)

      self.Move_Points(fit)
#      if loop > 1:
#        for det in chi_sq[loop]:
#          print det, "_________________________________"
#          for stat in chi_sq[loop][det]:
#            print "%.4f" % chi_sq[loop-1][det][stat] , "    ", \
#                  "%.4f" % chi_sq[loop][det][stat]
                  
      
    self.o_res.Write()
    self.o_prof.Write()

#########################################################################################
  # 
  def Draw_Line(self, seeds, cut):
    temp = {}
    for dim in ["x", "y"]:
      pos = dim+"_pos"
      a1, a2, a3, a4, c1, c2 = 0, 0, 0, 0 ,0 ,0
      for station in seeds:
        try:
          c1 += seeds[station][pos]*seeds[station]["z_pos"]
          c2 += seeds[station][pos]
          a1 += seeds[station]["z_pos"]**2
          a2 += seeds[station]["z_pos"]
          a4 += 1
        except KeyError:
          print "Fit error"
          print "Station ", station
          print seeds

      a = np.matrix([[a1,a2],[a2,a4]])
      b = np.matrix([[c1],[c2]])
      x = np.linalg.solve(a,b)
      m_dim = "m"+dim
      b_dim = "b"+dim
      temp[m_dim] = x[0]
      temp[b_dim] = x[1]
    return temp

#########################################################################################
  # 
  def Find_Coefficents(self, spaces, x_exp, y_exp, station):
    b11,b12,b13,b14 = 0,0,0,0
    b21,b22,b23,b24 = 0,0,0,0
    b31,b32,b33,b34 = 0,0,0,0
    b41,b42,b43,b44 = 0,0,0,0
    
    c1,c2,c3,c4 = 0,0,0,0
    d1,d2,d3,d4 = 0,0,0,0
    e1,e2,e3,e4 = 0,0,0,0
    
    temp = {}
    for event in spaces:
      seed = spaces[event][station]
      
      b11 += seed["x_pos"]**2
      b12 += seed["x_pos"] * seed["y_pos"]
      b13 += seed["x_pos"] * seed["z_pos"]
      b14 += seed["x_pos"]
      b22 += seed["y_pos"]**2
      b23 += seed["y_pos"] * seed["z_pos"]
      b24 += seed["y_pos"]
      b33 += seed["z_pos"]**2
      b34 += seed["z_pos"]
      b44 += 1

      c1 += seed["x_pos"] * x_exp[event]
      c2 += seed["y_pos"] * x_exp[event]
      c3 += seed["z_pos"] * x_exp[event]
      c4 += x_exp[event]

      d1 += seed["x_pos"] * y_exp[event]
      d2 += seed["y_pos"] * y_exp[event]
      d3 += seed["z_pos"] * y_exp[event]
      d4 += y_exp[event]

      e1 += seed["x_pos"] * seed["z_pos"]
      e2 += seed["y_pos"] * seed["z_pos"]
      e3 += seed["z_pos"] * seed["z_pos"]
      e4 += seed["z_pos"]

    b    = np.matrix([[b11,b12,b13,b14],[b12,b22,b23,b24],\
                      [b13,b23,b33,b34],[b14,b24,b34,b44]])
    c    = np.matrix([[c1],[c2],[c3],[c4]])
    d    = np.matrix([[d1],[d2],[d3],[d4]])
    e    = np.matrix([[e1],[e2],[e3],[e4]])
    row1 = np.linalg.solve(b,c)
    row2 = np.linalg.solve(b,d)
    row3 = np.linalg.solve(b,e)

    temp = np.matrix([[row1.item(0),row1.item(1),row1.item(2),row1.item(3)], \
                      [row2.item(0),row2.item(1),row2.item(2),row2.item(3)], \
                      [row3.item(0),row3.item(1),row3.item(2),row3.item(3)], \
                      [0,           0,           0,           1]])

    return temp

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
        
#        print "Old: ", self.position[detector][station]["x"]
#        print "Initial: ",x
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
#        print "New: ", self.position[detector][station]["x"]
#        print fit[detector][station],"\n"
        
#        step_fit = np.matrix([[m.cos(theta)*m.cos(phi), m.cos(theta)*m.sin(phi), \
#                              -m.sin(theta), x], \
#                              [m.sin(psi)*m.sin(theta)*m.cos(phi) - m.cos(psi)*m.sin(phi), \
#                               m.sin(psi)*m.sin(theta)*m.sin(phi) + m.cos(psi)*m.cos(phi), \
#                               m.cos(theta)*m.sin(psi), y], \
#                              [m.cos(psi)*m.sin(theta)*m.cos(phi) + m.sin(psi)*m.sin(phi), \
#                               m.cos(psi)*m.sin(theta)*m.sin(phi) - m.sin(psi)*m.cos(phi), \
#                               m.cos(theta)*m.cos(psi), z], \
#                              [0, 0, 0, 1]])        
        step_fit = np.matrix([[m.cos(sth)*m.cos(sph), m.cos(sth)*m.sin(sph), \
                              -m.sin(sth), sx], \
                              [m.sin(sps)*m.sin(sth)*m.cos(sph) - m.cos(sps)*m.sin(sph), \
                               m.sin(sps)*m.sin(sth)*m.sin(sph) + m.cos(sps)*m.cos(sph), \
                               m.cos(sth)*m.sin(sps), sy], \
                              [m.cos(sps)*m.sin(sth)*m.cos(sph) + m.sin(sps)*m.sin(sph), \
                               m.cos(sps)*m.sin(sth)*m.sin(sph) - m.sin(sps)*m.cos(sph), \
                               m.cos(sth)*m.cos(sps), sz], \
                              [0, 0, 0, 1]])

#        step_fit = np.matrix([[1,phi,-theta,x],[-phi+psi*theta,1+psi*phi*theta,psi,y],[theta+psi*phi,-psi+theta*phi,1,z],[0,0,0,1]])
#        print step_fit, "\n\n"

        for event in self.transit[detector]:
          point = self.transit[detector][event][station]
          tran = np.matrix([[point["x_pos"]],[point["y_pos"]], \
                            [point["z_pos"]],[1]])
          temp = step_fit*tran
#          temp = fit[detector][station]*tran

          point["x_pos"] += step * (temp.item(0) - point["x_pos"])
          point["y_pos"] += step * (temp.item(1) - point["y_pos"])
          point["z_pos"] += step * (temp.item(2) - point["z_pos"])
