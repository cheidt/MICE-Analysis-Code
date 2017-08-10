import math
from config_maus_analysis import cut_config as _config
import output_maus_analysis as _output
import numpy as np
import scipy
import scipy.stats
from scipy.stats import chisquare

class Cuts(object):
  def __init__ (self):
    self.o_time  = _output.Output("Time")
    self.o_emr   = _output.Output("EMR")
    self.o_other = _output.Output("Other")
    self.o_mass  = _output.Output("Mass")
    self.count   = {"Pass":0, "Total":0, "UMTot": 0, "UPTot": 0, "UETot": 0, "UUTot": 0, \
                    "MPass": 0, "PPass": 0, "EPass": 0, "UPass": 0, "TEMR": 0, "MMSDR": 0, "PMSDR": 0, \
                    "EMSDR": 0, "UMSDR": 0, "DMTot": 0, "DPTot": 0, "DETot": 0, "DUTot": 0}

  def Write(self):
    # self.o_time.Gaus_Fit(27.0, 29.0, "DTime")
    # self.o_time.Gaus_Fit(29.0, 31.5, "DTime")
    # self.o_time.Gaus_Fit(31.5, 34.0, "DTime")
    #self.o_emr.Add_Graphs("Chi2_Graph", "EMR Chi2", "Graph_EMR_Entire_XY_Ch2_muon", \
    #                       "Graph_EMR_Entire_XY_Ch2_pion", "Graph_EMR_Entire_XY_Ch2_electron")
    self.o_time.Write()
    self.o_emr.Write()
    self.o_other.Write()
    self.o_mass.Write()

  def Out_Count(self):
    return self.count

###############################################################################
# Runs particle cuts on an event by event basis.  These should be single
#  particle events.  Cuts include checks for single a single particle track in
#  the TOFs and Trackers, fidelity in the Trackers, and PID in TOFs and the EMR
###############################################################################
  def Process_Cuts(self, data):
    self.count["Total"] += 1
    pass_cut = {"Rad_Diff": [False, [-10], [-10]], "UTOF_Time": [False, [-10], -10], "Only_One_TOF0": False, \
                "Only_One_TOF1": False, "Mass_Cut": [False, [-10]], "Only_One_UTrack": False, \
                "Only_One_DTrack": False, "UTracks": [False, [-10]], "TOF0": [False, 0, [-10]], \
                "TOF1": [False, 0, [-10]], "TOF2": [False, 0, [-10]], "Only_One_TOF2": False, \
                "Only_One_UTOF": False, "Only_One_UTOF_Track": False, "Only_One_All": False, \
                "Only_One_Track": False, "DTOF_Time": [False, [-10]], "DTracks": [False, -10], \
                "Type": ["unknown", "unknown", "unknown"], "UTrack_Good": [False, [False]], \
                "DTrack_Good": [False, [False]], "Only_One_DTOF": False, "EMR_Chi2": [False, -10, -10, -10], \
                "EMR_Den": [False, -10], "UTracker_Scraping": [False, -10], "DTracker_Scraping": [False, -10], \
                "Pass": False, "Mom_Cut": [[False], [False]]}

    pass_cut = self.Check_TOF(data, pass_cut)
    pass_cut = self.Check_Tracker(data, pass_cut)
    pass_cut = self.Check_EMR(data, pass_cut)
    pass_cut = self.Logic_Tests(pass_cut)
    pass_cut = self.Check_Mass(data, pass_cut)
    if pass_cut["UTracks"][0] and pass_cut["TOF0"][0] and pass_cut["TOF1"][0]:
      pass_cut["Mass_Cut"] = self.Mass_Cut(pass_cut["TOF0"][2], pass_cut["TOF1"][2], data["tracker_tracks"]["upstream"], len(data["tracker_tracks"]["upstream"]))
    self.Cut_Prints(pass_cut)
#    raw_input("Press Enter to Exit")
    return pass_cut


  def Check_TOF(self, data, pass_cut):
    if "TOF0_space_points" in data:
      pass_cut["TOF0"] = [True, len(data["TOF0_space_points"]), [-10]]
      if pass_cut["TOF0"][1] == 1:
        pass_cut["Only_One_TOF0"] = [True]
      TOF0Time = []
      for i in range(pass_cut["TOF0"][1]):
        TOF0Time.append(data["TOF0_space_points"][i]["time"])
      pass_cut["TOF0"][2] = TOF0Time

    if "TOF1_space_points" in data:
      pass_cut["TOF1"] = [True, len(data["TOF1_space_points"]), [-10]]
      if pass_cut["TOF1"][1] == 1:
        pass_cut["Only_One_TOF1"] = True
      TOF1Time = []
      for j in range(pass_cut["TOF1"][1]):
        TOF1Time.append(data["TOF1_space_points"][j]["time"])
      pass_cut["TOF1"][2] = TOF1Time

    if "TOF2_space_points" in data:
      pass_cut["TOF2"] = [True, len(data["TOF2_space_points"]), [-10]]
      if pass_cut["TOF2"][1] == 1:
        pass_cut["Only_One_TOF2"] = True
      TOF2Time = []
      for k in range(pass_cut["TOF2"][1]):
        TOF2Time.append(data["TOF2_space_points"][k]["time"])
      pass_cut["TOF2"][2] = TOF2Time

    if pass_cut["Only_One_TOF1"] and pass_cut["Only_One_TOF0"]:
      pass_cut["Only_One_UTOF"] = True

    if pass_cut["Only_One_TOF1"] and pass_cut["Only_One_TOF2"]:
      pass_cut["Only_One_DTOF"] = True

    if pass_cut["TOF1"][0] and pass_cut["TOF2"][0]:
      pass_cut["DTOF_Time"] = self.TOF_Timing_Cut(TOF1Time, TOF2Time)
      tof_time = pass_cut["DTOF_Time"][1][0]
      # print tof_time
      mu = _config["D_c_mu"] * math.e ** (-0.5 * ((tof_time - _config["D_mu_mu"]) / _config["D_sigma_mu"]) ** 2)
      pi = _config["D_c_pi"] * math.e ** (-0.5 * ((tof_time - _config["D_mu_pi"]) / _config["D_sigma_pi"]) ** 2)
      el = _config["D_c_el"] * math.e ** (-0.5 * ((tof_time - _config["D_mu_el"]) / _config["D_sigma_el"]) ** 2)
      if mu + pi + el > 0.0001:
        pass_cut["DTOF_Time"][2] = mu / (mu + pi + el)
        Dpion_ratio = pi / (mu + pi + el)
        Delec_ratio = el / (mu + pi + el)
        if pass_cut["DTOF_Time"][2] > _config["particle_selection_ratio"]:
          pass_cut["Type"][0] = "muon"
          self.count["DMTot"] += 1
        elif Dpion_ratio > _config["particle_selection_ratio"]:
          pass_cut["Type"][0] = "pion"
          self.count["DPTot"] += 1
        elif Delec_ratio > _config["particle_selection_ratio"]:
          pass_cut["Type"][0] = "electron"
          self.count["DETot"] += 1
        else:
          self.count["DUTot"] += 1

    if pass_cut["TOF0"][0] and pass_cut["TOF1"][0]:
      pass_cut["UTOF_Time"] = self.TOF_Timing_Cut(TOF0Time, TOF1Time)
      tof_time = pass_cut["UTOF_Time"][1][0]
      mu = _config["U_c_mu"] * math.e**(-0.5 * ((tof_time - _config["U_mu_mu"]) / _config["U_sigma_mu"])**2)
      pi = _config["U_c_pi"] * math.e**(-0.5 * ((tof_time - _config["U_mu_pi"]) / _config["U_sigma_pi"])**2)
      el = _config["U_c_el"] * math.e**(-0.5 * ((tof_time - _config["U_mu_el"]) / _config["U_sigma_el"])**2)
      # print mu, pi, el
      if mu + pi + el > 0.00001:
        Umuon_ratio = mu / (mu + pi + el)
        Upion_ratio = pi / (mu + pi + el)
        Uelec_ratio = el / (mu + pi + el)
        if Umuon_ratio > _config["particle_selection_ratio"]:
          pass_cut["Type"][0] = "muon"
          pass_cut["UTOF_Time"][2] = Umuon_ratio
          self.count["UMTot"] += 1
        elif Upion_ratio > _config["particle_selection_ratio"]:
          pass_cut["Type"][0] = "pion"
          pass_cut["UTOF_Time"][2] = Upion_ratio
          self.count["UPTot"] += 1
        elif Uelec_ratio > _config["particle_selection_ratio"]:
          pass_cut["Type"][0] = "electron"
          pass_cut["UTOF_Time"][2] = Uelec_ratio
          self.count["UETot"] += 1
        else:
          self.count["UUTot"] += 1
    return pass_cut


  def Check_Tracker(self, data, pass_cut):
    if "tracker_tracks" in data:
      if len(data["tracker_tracks"]["upstream"]) > 0:
        good = []; mom_cut =[]
        pass_cut["UTracks"] = [True, len(data["tracker_tracks"]["upstream"])]
        if pass_cut["UTracks"][1] == 1:
          pass_cut["Only_One_UTrack"] = True
        for i in range(pass_cut["UTracks"][1]):
          good.append(data["tracker_tracks"]["upstream"][i]["good"])
          if good[-1]:
            pass_cut["UTrack_Good"][0] = True
          for point in data["tracker_tracks"]["upstream"][i]["track_points"]:
            if point["station"] == 5 and point["plane"]  == 0:
              if abs(point["x_mom"]) > .0001 and abs(point["y_mom"]) > .0001:
                mom_cut.append(True)
              else:
                mom_cut.append(False)
        pass_cut["Mom_Cut"][0] = mom_cut
        pass_cut["UTrack_Good"][1] = good

      if len(data["tracker_tracks"]["downstream"]) > 0:
        good = []; mom_cut = []
        pass_cut["DTracks"] = [True, len(data["tracker_tracks"]["downstream"])]
        if pass_cut["DTracks"][1] == 1:
          pass_cut["Only_One_DTrack"] = True
        for i in range(pass_cut["DTracks"][1]):
          good.append(data["tracker_tracks"]["downstream"][i]["good"])
          if good[-1]:
            pass_cut["DTrack_Good"][0] = True
          for point in data["tracker_tracks"]["downstream"][i]["track_points"]:
            if point["station"] == 5 and point["plane"]  == 0:
              if abs(point["x_mom"]) > .0001 and abs(point["y_mom"]) > .0001:
                mom_cut.append(True)
              else:
                mom_cut.append(False)
        pass_cut["Mom_Cut"][1] = mom_cut
        pass_cut["DTrack_Good"][1] = good

    if pass_cut["UTracks"][0]:
      pass_cut["Rad_Diff"] = self.Radial_Diffuser_Cut(data["tracker_tracks"]["upstream"], \
                                                        pass_cut["UTracks"][1])
      pass_cut["UTracker_Scraping"] = self.Test_Tracker_Scraping(data["tracker_tracks"]["upstream"])

    if pass_cut["DTracks"][0]:
      pass_cut["DTracker_Scraping"] = self.Test_Tracker_Scraping(data["tracker_tracks"]["downstream"])
    return pass_cut


  def Check_EMR(self, data, pass_cut):
    # name = "EMR_SP_XZ_MC"
    # self.o_emr.Delete(name)
    # name = "EMR_SP_XZ_SC"
    # self.o_emr.Delete(name)
    # name = "EMR_SP_YZ_MC"
    # self.o_emr.Delete(name)
    # name = "EMR_SP_YZ_SC"
    # self.o_emr.Delete(name)
    if "EMR_tracks" in data:
      self.count["TEMR"] += 1
      data = data["EMR_tracks"]
      # pass_cut["EMR_Chi2"] = self.EMR_Chi2_Test(data)
      # pass_cut["EMR_Den"]  = self.EMR_Density_Test(data)
      for track in data:
        if track["type"] == "mother":
          pass_cut["EMR_Den"][1]  = track["mden_ratio"]
          pass_cut["EMR_Chi2"][3] = track["chi2"]
          if track["mden_ratio"] > .9:
            pass_cut["EMR_Den"][0] = True
          if track["chi2"] < 5:
            pass_cut["EMR_Chi2"][0] = True
      pass_cut["EMR_Mass"] = self.EMR_Mass_Test(data, pass_cut)
      # for track in data:
      #   name = "EMR_Time_" + pass_cut["Type"][0]
      #   title = "EMR Global Timing " + pass_cut["Type"][0]
      #   self.o_emr.Fill(name, title, track["time"], 300, 0, 3000000)
      #
      #   name = "EMR_MCharge_Density_" + pass_cut["Type"][0]
      #   title = "EMR Mother Charge Density " + pass_cut["Type"][0]
      #   self.o_emr.Fill(name, title, track["mcharge_ratio"], 200, 0, 5)
      #
      #   name = "EMR_SDensity_" + pass_cut["Type"][0]
      #   title = "EMR Secondary Plane Density " + pass_cut["Type"][0]
      #   self.o_emr.Fill(name, title, track["sden_ratio"], 200, 0, 1.2)
      #
      #   name = "EMR_Range_" + pass_cut["Type"][0]
      #   title = "EMR Range " + pass_cut["Type"][0]
      #   self.o_emr.Fill(name, title, track["range"], 100, 0, 900)
      #
      #   name = "EMR_Momentum_" + pass_cut["Type"][0]
      #   title = "EMR Momentum " + pass_cut["Type"][0]
      #   self.o_emr.Fill(name, title, track["mom"], 300, 5, 300)
    return(pass_cut)


  def Logic_Tests(self, pass_cut):
    if pass_cut["Only_One_UTrack"] and pass_cut["Only_One_DTrack"] and \
       pass_cut["UTrack_Good"][0] and pass_cut["DTrack_Good"][0]:
      pass_cut["Only_One_Track"] = True

    if pass_cut["Only_One_Track"] and pass_cut["Only_One_UTOF"]:
      pass_cut["Only_One_UTOF_Track"] = True

    if pass_cut["Only_One_UTOF_Track"] and pass_cut["Only_One_TOF2"]:
      pass_cut["Only_One_All"] = True

    if pass_cut["Only_One_All"] and pass_cut["UTOF_Time"][0] and pass_cut["UTracker_Scraping"][0] and \
       pass_cut["DTracker_Scraping"][0] and pass_cut["Mom_Cut"][0][0] and pass_cut["Mom_Cut"][1][0]:
      #and pass_cut["Mass_Cut"][0]:
      self.count["Pass"] += 1
      if pass_cut["Type"][0] == "muon":
        self.count["MPass"] += 1
        pass_cut["Pass"] = True
      if pass_cut["Type"][0] == "pion":
        self.count["PPass"] += 1
      if pass_cut["Type"][0] == "electron":
        self.count["EPass"] += 1
      if pass_cut["Type"][0] == "unknown":
        self.count["UPass"] += 1

    return pass_cut


  def Check_Mass(self, data, pass_cut):
    me = _config["me"]; mm = _config["mm"]; mp = _config["mp"]; density = _config["density"]
    an = _config["atomic_number"]; am = _config["atomic_mass"]; K = _config["K"]
    I  = _config["mean_excitation"] * an
    dz = 4.80
    if pass_cut["Only_One_All"] and not pass_cut["Type"][0] == "unknown":
      for detector in ["upstream", "downstream"]:
        key = str(detector[0].capitalize() + "TOF_Time")
        ptime = pass_cut[key][1][0]
        ctime = _config["{}_timing".format(detector)]
        beta  = ctime/ptime

        if beta > 0 and beta < 1:
          for point in data["tracker_tracks"][detector][0]["track_points"]:
            if point["plane"] == 0 and point["station"] == 5:
              momentum = point
          gamma = math.sqrt(1. - beta ** 2.)
          mpart = mm
          T_max = 2. * me * beta**2. * gamma**2. / (1.+(2.*gamma * me / mpart) + (me / mpart)**2.)
          dedl = K * (an / am) / beta**2. * (0.5*math.log(2.*me * beta**2. * gamma**2. * T_max / I**2.) - beta**2)
          de = dedl * dz * 100 * density

          if de > 0:
            loss_mom = math.sqrt((de + mpart)**2. - mpart**2.)
            restored_mom = math.sqrt(momentum["z_mom"]**2. + momentum["x_mom"]**2. + momentum["y_mom"]**2.) + loss_mom
            raw_mass = restored_mom / beta
            final_mass = gamma * raw_mass

            mu = _config["{}m_c_mu".format(detector[0])] * math.e**(-0.5*((final_mass - \
                 _config["{}m_mu_mu".format(detector[0])]) / _config["{}m_sigma_mu".format(detector[0])])**2)
            pi = _config["{}m_c_pi".format(detector[0])] * math.e ** (-0.5 * ((final_mass - \
                 _config["{}m_mu_pi"].format(detector[0])) / _config["{}m_sigma_pi"].format(detector[0]))**2)
            el = _config["{}m_c_el".format(detector[0])] * math.e**(-0.5*((final_mass - \
                 _config["{}m_mu_el"].format(detector[0])) / _config["{}m_sigma_el"].format(detector[0]))**2)
            # print mu, pi, el
            if mu + pi + el > 0.00001:
              muon_ratio = mu / (mu + pi + el)
              pion_ratio = pi / (mu + pi + el)
              elec_ratio = el / (mu + pi + el)

              if detector == "upstream":
                i = 1
              else:
                i = 2
              if muon_ratio > _config["particle_selection_ratio"]:
                pass_cut["Type"][i] = "muon"
              if pion_ratio > _config["particle_selection_ratio"]:
                pass_cut["Type"][i] = "pion"
              if elec_ratio > _config["particle_selection_ratio"]:
                pass_cut["Type"][i] = "electron"

            name  = "Mass_{}".format(detector)
            title = "Mass {}".format(detector)
            self.o_mass.Fill(name, title, final_mass, 600, 0, 300)
    return pass_cut


  def Mass_Cut(self, time0, time1, tracks, size):
    #dz = 7.64231643643#dz = 7.94231643643
    dz = 7.610
    t_e = 24.75
    m_e = 0.510999
    m_m = 105.66
    final_mass = []
    for t0 in time0:
      for t1 in time1:
        for i in range(size):
          for l in range(len(tracks[i]["track_points"])):
            if tracks[i]["track_points"][l]["station"] == 5 and \
                            tracks[i]["track_points"][l]["plane"] == 0:
              track = tracks[i]["track_points"][l]
          time = (t1 - t0)
          beta = t_e/time
          density = .001204
          I = 11.6 * 7.311 / 1000000

          #new_mass = False
          if beta < 1 and beta > 0:
            gamma = math.sqrt(1.0 - beta**2)
            T_max = 2*m_e*beta**2*gamma**2/(1+(2*gamma*m_e/m_m)+(m_e/m_m)**2)
            energy_loss_per_step = 0.307075 * (7.311/14.666) / beta**2 * (0.5 * math.log(2.*m_e*beta**2*gamma**2*T_max/I**2) - beta**2)
            energy_loss = energy_loss_per_step * dz * 100 * density

            if energy_loss > 0:
              loss_mom = ((energy_loss + m_m)**2 - m_m**2)**(.5)
              restored_mom = math.sqrt(track["z_mom"]**2 + track["x_mom"]**2 + track["y_mom"]**2) + loss_mom
              raw_mass = restored_mom / beta
              final_mass.append(gamma * raw_mass)
    for mass in final_mass:
      if abs(mass - 103.22) < 7.8:
        return [True, final_mass]
    return [False, final_mass]

            #else:
            #new_mass = True

            #if new_mass == True and energy_loss < 0:
            #T_max = 2*m_e*beta**2*gamma**2/(1+(2*gamma*m_e/m_p)+(m_e/m_p)**2)
            #energy_loss_per_step = 0.307075 * (7.311/14.666) / beta**2 * (0.5 * math.log(2.*m_e*beta**2*gamma**2*T_max/I**2) - beta**2)
            #energy_loss = energy_loss_per_step * dz * 100 * density

            #if energy_loss > 0:
            #loss_mom = ((energy_loss + m_p)**2 - m_p**2)**(.5)
            #restored_mom = math.sqrt(track["z_mom"]**2 + track["x_mom"]**2 + track["y_mom"]**2) + loss_mom
            #raw_mass = restored_mom / beta
            #final_mass = gamma * raw_mass

            #if final_mass > 500.:
            #name = "Mass_Plot"
            #title = "Mass"
            #self.o_time.Fill(name, title, final_mass, 1200, 0, 1200)
            #print final_mass
            #else:
            #print "Whoops"
            #return False


  def TOF_Timing_Cut(self, tof1, tof2):
    time = []
    for time1 in tof1:
      for time2 in tof2:
        time.append(time2 - time1)
    for dif in time:
      if dif < 100 and dif > 20:
        return [True, time, -10]
    return [False, time, -10]


  def Radial_Diffuser_Cut(self, tracks, size):
    radius = []
    z_p    = []
    for i in range(size):
      for j in range(len(tracks[i]["track_points"])):
        track_pt = tracks[i]["track_points"][j]
        if track_pt["station"] == 5 and track_pt["plane"] == 0:
          x_i = track_pt["x_pos"]
          y_i = track_pt["y_pos"]
          x_p = track_pt["x_mom"]
          y_p = track_pt["y_mom"]
          z_p.append(track_pt["z_mom"])
          d_z = 800 #mm

          x_f = (d_z*x_p/z_p[-1])+x_i
          y_f = (d_z*y_p/z_p[-1])+y_i
          radius.append(math.sqrt(x_f**2 + y_f**2))

    for rad in radius:
      if rad < 90:
        return [True, z_p, radius]
    return [False, z_p, radius]


  def Test_Tracker_Scraping(self, data):
    length = -10
    for track in data:
      if len(track["prtrks"]) > 1:
        print "Should only be one track:", len(track["prtrks"])
        raw_input("Press Enter to Exit")
      prtrack = track["prtrks"][0]
      radius  = prtrack["radius"]
      xo      = prtrack["circle_x0"]
      yo      = prtrack["circle_y0"]
      length  = math.sqrt(xo**2 + yo**2) + radius
      if length < 150.0:
        return [True, length]
    return [False, length]


  ###############################################################################
  # Pulls out EMR space point info from both mother and secondary tracks and
  # preforms a line fit and chi2 analysis.
  ###############################################################################
  def EMR_Chi2_Test(self, data):
    xz_array = [[],[]]; yz_array = [[],[]]
    x_chi = -10; y_chi = -10; tot_chi = -10
    for track in data:
      if len(track["space_points"]) > 0:
        for sp in range(len(track["space_points"])):
          space = track["space_points"][sp]
          if not space["x_pos"] in xz_array[0] and not space["z_pos"] in xz_array[1]:
            xz_array[0].append(space["x_pos"])
            xz_array[1].append(space["z_pos"])
          if not space["y_pos"] in yz_array[0] and not space["z_pos"] in yz_array[1]:
            yz_array[0].append(space["y_pos"])
            yz_array[1].append(space["z_pos"])

    if len(xz_array[1]) > 5 and len(yz_array[1]) > 5:
      x_expect = []; y_expect = []
      x_fit = np.polyfit(xz_array[1], xz_array[0], 1, full=True)
      y_fit = np.polyfit(yz_array[1], yz_array[0], 1, full=True)
      for i in range(len(xz_array[1])):
        x_expect.append(xz_array[1][i] * x_fit[0][0] + x_fit[0][1])
      for i in range(len(yz_array[1])):
        y_expect.append(yz_array[1][i] * y_fit[0][0] + y_fit[0][1])
        # print "Expected Y Values"
        # print "y = ", y_fit[0][0], "*", yz_array[1][i], " + ", y_fit[0][1]
        # print yz_array[0][i], " / ", y_expect[i], "\n\n"
      x_chi = chisquare(xz_array[0], x_expect)
      y_chi = chisquare(yz_array[0], y_expect)
      tot_chi = math.sqrt(x_chi[0] ** 2 + y_chi[0] ** 2)
      if tot_chi < 10:
        return [True, abs(x_chi[0]), abs(y_chi[0]), tot_chi]
      else:
        return [False, abs(x_chi[0]), abs(y_chi[0]), tot_chi]
    return [False, x_chi, y_chi, tot_chi]

  ###############################################################################
  # Pulls EMR plane hits.  Determines the plane furthest back hit and the
  # total number of plane hits to find density of planes activated in EMR
  ###############################################################################
  def EMR_Density_Test(self, data):
    plane = []
    density = -10
    for track in data:
      if len(track["plane_hits"]) > 0:
        for ph in range(len(track["plane_hits"])):
          plane.append(track["plane_hits"][ph]["plane"])
    if len(plane) > 2:
      plane = list(set(plane))
      final = float(plane[-1]+1.0)
      density = float(len(plane))/final
      if density > .90:
        return [True, density]
    return[False, density]


  def EMR_Mass_Test(self, data, pass_cut):
    q1 = 0.007297**(0.5)
    c = .001239842/1000000.0
    # if range(len(data)) > 1:
      # print "Two EMR tracks"
    # print pass_cut["Type"]
    mass_ma = data[0]["charge_MA"] / q1 * (data[0]["charge_MA"] / (q1 * data[0]["range"]) - 1)
    mass_sa = data[0]["charge_SA"] / q1 * (data[0]["charge_SA"] / (q1 * data[0]["range"]) - 1)
    #
    # print "Charge MA: ", data[0]["charge_MA"]
    # print "Mass MA:   ", mass_ma
    # print "Charge SA: ", data[0]["charge_SA"]
    # print "Mass SA:   ", mass_sa
    # print "Range:     ", data[0]["range"], "\n"

    mass_ma = data[0]["charge_MA"] / q1 * (data[0]["charge_MA"] / (q1 * data[0]["range"] * c) - 1)
    mass_sa = data[0]["charge_SA"] / q1 * (data[0]["charge_SA"] / (q1 * data[0]["range"] * c) - 1)
    #
    # print "After Range Conversion: "
    # print "Mass MA:   ", mass_ma
    # print "Mass SA:   ", mass_sa, "\n"

    q2 = 1.60217662 * 10.0**(-19)
    k  = 8.987552   * 10.0**(9)
    c  = 5.60958    * 10.0**(29)
    mass_ma = data[0]["charge_MA"] / q2 * (data[0]["charge_MA"] / (q2 * k * data[0]["range"] / 1000) - 1)
    mass_sa = data[0]["charge_SA"] / q2 * (data[0]["charge_SA"] / (q2 * k * data[0]["range"] / 1000) - 1)
    #
    # print "Using MKS Units: "
    # print "Mass MA:   ", mass_ma
    # print "Mass SA:   ", mass_sa, "\n\n\n"

  ###############################################################################
  # Prints out everything I care to look at
  ###############################################################################
  def Cut_Prints(self, data):
    name = "Cut_Counts"
    title = "Number of events to pass cut"
    number = 0
    self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["TOF0"][0]:
      number = 1
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["TOF1"][0]:
      number = 2
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["TOF2"][0]:
      number = 3
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["Only_One_TOF0"]:
      number = 4
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["Only_One_TOF1"]:
      number = 5
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["Only_One_TOF2"]:
      number = 6
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["UTOF_Time"]:
      number = 7
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["Mass_Cut"]:
      number = 8
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["UTracks"]:
      number = 9
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["DTracks"]:
      number = 10
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["Only_One_UTrack"]:
      number = 11
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["Only_One_DTrack"]:
      number = 12
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["UTrack_Good"]:
      number = 13
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["DTrack_Good"]:
      number = 14
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["UTracker_Scraping"]:
      number = 15
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["DTracker_Scraping"]:
      number = 16
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["Only_One_UTOF"]:
      number = 17
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["Only_One_Track"]:
      number = 18
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["Only_One_UTOF_Track"]:
      number = 19
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["Only_One_All"]:
      number = 20
      self.o_other.Fill(name, title, number, 22, 0, 22)
    if data["Pass"]:
      number = 21
      self.o_other.Fill(name, title, number, 22, 0, 22)

    for i in range(len(data["Mass_Cut"][1])):
      name = "Mass_Plot"
      title = "Mass"
      self.o_mass.Fill(name, title, data["Mass_Cut"][1][i], 500, 0, 500)

    for i in range(len(data["DTOF_Time"][1])):
      name = "DTime"
      title = "Downstream Timing"
      self.o_time.Fill(name, title, data["DTOF_Time"][1][i], 500, 25, 45)

      name = "DTime_" + data["Type"][0]
      title = "Downstream Timing " + data["Type"][0]
      self.o_time.Fill(name, title, data["DTOF_Time"][1][i], 500, 25, 45)


    for i in range(len(data["UTOF_Time"][1])):
      name = "UTime"
      title = "Upstream Timing"
      self.o_time.Fill(name, title, data["UTOF_Time"][1][i], 500, 20, 40)

      if data["UTOF_Time"][2] > 0.0:
        ratio = data["UTOF_Time"][2]
        name = "UTime_Weighted_" + data["Type"][0]
        title = "Upstream Timing Weighted for " + data["Type"][0]
        self.o_time.Fill(name, title, data["UTOF_Time"][1][i], 500, 20, 40, weight = ratio)
      else:
        name = "UTime_" + data["Type"][0]
        title = "Upstream Timing " + data["Type"][0]
        self.o_time.Fill(name, title, data["UTOF_Time"][1][i], 500, 20, 40)
      if not data["Type"][0] == "unknown":
        name = "UTime_" + data["Type"][0]
        title = "Upstream Timing " + data["Type"][0]
        self.o_time.Fill(name, title, data["UTOF_Time"][1][0], 500, 20 , 40)

    if data["Only_One_UTOF_Track"]:
      name = "Diff_Cut"
      title = "TOF Timing vs Track Momentum"
      if len(data["Rad_Diff"][1]) > 1:
        raw_input("Press Enter to Exit")
      time = data["UTOF_Time"][1][0]
      mom  = data["Rad_Diff"][1][0]
      self.o_time.Fill(name, title, mom, time, \
                      250, 0, 400, 250, 20, 60)

    if data["Only_One_UTOF_Track"] and data["Rad_Diff"]:
      name = "Diff_Cut_True"
      title = "TOF Timing vs Track Momentum After Cut"
      time = data["UTOF_Time"][1][0]
      mom  = data["Rad_Diff"][1][0]
      self.o_time.Fill(name, title, mom, time, \
                      250, 0, 400, 250, 20, 60)

    if data["EMR_Den"][1] > -1:
      name = "EMR_MDensity_" + data["Type"][0]
      title = "EMR Plane Density " + data["Type"][0]
      self.o_emr.Fill(name, title, data["EMR_Den"][1], 200, 0, 1.2)

      name = "EMR_Ch2_" + data["Type"][0]
      title = "EMR Track Chi2" + data["Type"][0]
      self.o_emr.Fill(name, title, data["EMR_Chi2"][3], 100, 0, 200)

    # if data["EMR_Chi2"][2] > -1:
    #   name = "EMR_Entire_X_Ch2_" + data["Type"][0]
    #   title = "EMR All Tracks X Chi2 " + data["Type"][0]
    #   self.o_emr.Fill(name, title, data["EMR_Chi2"][1], 250, 0, 20000)
    #
    #   name = "EMR_Entire_Y_Ch2_" + data["Type"][0]
    #   title = "EMR All Tracks Y Chi2 " + data["Type"][0]
    #   self.o_emr.Fill(name, title, data["EMR_Chi2"][2], 250, 0, 20000)
    #
    #   name = "Narrow_EMR_Entire_X_Ch2_" + data["Type"][0]
    #   title = "EMR All Tracks X Chi2 " + data["Type"][0]
    #   self.o_emr.Fill(name, title, data["EMR_Chi2"][1], 250, 0, 200)
    #
    #   name = "Narrow_EMR_Entire_Y_Ch2_" + data["Type"][0]
    #   title = "EMR All Tracks Y Chi2 " + data["Type"][0]
    #   self.o_emr.Fill(name, title, data["EMR_Chi2"][2], 250, 0, 200)
    #
    #   name = "EMR_Entire_XY_Ch2_" + data["Type"][0]
    #   title = "EMR All Tracks XY Chi2 " + data["Type"][0]
    #   self.o_emr.Fill(name, title, data["EMR_Chi2"][1], data["EMR_Chi2"][2], 250, 0, 20000, 500, 0, 20000)
    #
    #   name = "Graph_EMR_Entire_XY_Ch2_" + data["Type"][0]
    #   title = "EMR All Tracks XY Chi2 " + data["Type"][0]
    #   col = 1
    #   if data["Type"][0] == "muon":
    #     col = 2
    #   if data["Type"][0] == "pion":
    #     col = 4
    #   if data["Type"][0] == "electron":
    #     col = 8
    #   self.o_emr.Graph(name, title, data["EMR_Chi2"][1], data["EMR_Chi2"][2], color = col)
    #
    #   name = "Narrow_EMR_Entire_XY_Ch2_" + data["Type"][0]
    #   title = "EMR All Tracks XY Chi2 " + data["Type"][0]
    #   self.o_emr.Fill(name, title, data["EMR_Chi2"][1], data["EMR_Chi2"][2], 250, 0, 200, 500, 0, 200)
    #
    #   name = "EMR_Entire_T_Ch2_" + data["Type"][0]
    #   title = "EMR All Tracks T Chi2 " + data["Type"][0]
    #   self.o_emr.Fill(name, title, data["EMR_Chi2"][3], 250, 0, 20000)
    #
    #   name = "Narrow_EMR_Entire_T_Ch2_" + data["Type"][0]
    #   title = "EMR All Tracks T Chi2 " + data["Type"][0]
    #   self.o_emr.Fill(name, title, data["EMR_Chi2"][3], 250, 0, 200)
    #
    # if data["EMR_Den"][1] > -1:
    #   name = "EMR_Plane_Density_" + data["Type"][0]
    #   title = "EMR All Plane Density " + data["Type"][0]
    #   self.o_emr.Fill(name, title, data["EMR_Den"][1], 150, 0, 1.5)