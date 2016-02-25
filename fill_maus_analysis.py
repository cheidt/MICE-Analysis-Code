### Error dealling with vector<double> in at least track points
###  in error and covariance.  Output ex:
###  <ROOT.vector<double> object at 0x409f120>
###  Too tried to deal with it tonight.

import config_maus_analysis as _config

#########################################################################################
  # Control function called directly from process routines.  Reads in single
  #   MICE event in MAUS format and returns nestled dictionary with all event
  #   information.
def Fill_from_Data(recon, mc, is_recon, is_mc):
  results = {}
  if is_recon:
    if len(recon.GetSciFiEvent().digits())               != 0:
      results["tracker_digits"]       = Digits(recon.\
                                        GetSciFiEvent().digits())
    if len(recon.GetSciFiEvent().clusters())             != 0:
      results["tracker_clusters"]     = Clusters(recon.\
                                        GetSciFiEvent().clusters())
    if len(recon.GetSciFiEvent().spacepoints())          != 0:
      results["tracker_space_points"] = Space_Points(recon.\
                                        GetSciFiEvent().spacepoints())
    if len(recon.GetSciFiEvent().straightprtracks())     != 0:
      results["tracker_straight_pr"]  = Straight_Pattern_Recon(recon.\
                                        GetSciFiEvent().straightprtracks())
    if len(recon.GetSciFiEvent().helicalprtracks())      != 0:
      results["tracker_helical_pr"]   = Helical_Pattern_Recon(recon.\
                                        GetSciFiEvent().helicalprtracks())
    if len(recon.GetSciFiEvent().scifitracks())          != 0:
      results["tracker_tracks"]       = Tracks(recon.\
                                        GetSciFiEvent().scifitracks())
    if len(recon.GetTOFEvent().GetTOFEventSpacePoint().\
                               GetTOF0SpacePointArray()) != 0:
      results["TOF0_space_points"]    = TOF(recon.GetTOFEvent().\
                                                  GetTOFEventSpacePoint().\
                                                  GetTOF0SpacePointArray())
    if len(recon.GetTOFEvent().GetTOFEventSpacePoint().\
                               GetTOF1SpacePointArray()) != 0:
      results["TOF1_space_points"]    = TOF(recon.GetTOFEvent().\
                                                  GetTOFEventSpacePoint().\
                                                  GetTOF1SpacePointArray())
    if len(recon.GetTOFEvent().GetTOFEventSpacePoint().\
                               GetTOF2SpacePointArray()) != 0:
      results["TOF2_space_points"]    = TOF(recon.GetTOFEvent().\
                                                  GetTOFEventSpacePoint().\
                                                  GetTOF2SpacePointArray())
  if is_mc == True:
    if len(mc.GetVirtualHits()) != 0:
      results["virtual_points"] = Virtual_Hits(mc.GetVirtualHits())

  return results

#########################################################################################
  # Fills tracker digits
def Digits(digits):
  container = {"upstream":{1:[], 2:[], 3:[], 4:[], 5:[]}, \
             "downstream":{1:[], 2:[], 3:[], 4:[], 5:[]}}
  for di in range(len(digits)):
    digit    = fill_digits(digits[di])
    tracker  = digit["tracker"]
    station  = digit["station"]
    detector = "upstream" if tracker == 0 else "downstream"
    container[detector][station].append(digit)
  return container

#########################################################################################
  # Fills tracker clusters
def Clusters(clusters):
  container = {"upstream":{1:[], 2:[], 3:[], 4:[], 5:[]}, \
             "downstream":{1:[], 2:[], 3:[], 4:[], 5:[]}}
  for cl in range(len(clusters)):
    cluster  = fill_clusters(clusters[cl])
    tracker  = cluster["tracker"]
    station  = cluster["station"]
    detector = "upstream" if tracker == 0 else "downstream"
    container[detector][station].append(cluster)
  return container

#########################################################################################
  # Collects information from tracker space points. 
def Space_Points(spaces):
  container = {"upstream":{1:[], 2:[], 3:[], 4:[], 5:[]}, \
             "downstream":{1:[], 2:[], 3:[], 4:[], 5:[]}}
  for sp in range(len(spaces)):
    space    = fill_space_points(spaces[sp])
    tracker  = space["tracker"]
    station  = space["station"]
    detector = "upstream" if tracker == 0 else "downstream"
    container[detector][station].append(space)
  return container

#########################################################################################
  # Collects straight track information from pattern recognition.
def Straight_Pattern_Recon(prtrks):
  container = {"upstream":[], "downstream":[]}
  for st in range(len(prtrks)):
    st_prtr  = fill_straight_pattern_recon(prtrks[st])
    tracker  = st_prtr["tracker"]
    detector = "upstream" if tracker == 0 else "downstream"
    container[detector].append(st_prtr)
  return container

#########################################################################################
  # Collects helical track information from pattern recognition.
def Helical_Pattern_Recon(prtrks):
  container = {"upstream":[], "downstream":[]}
  for ht in range(len(prtrks)):
    hl_prtr = fill_helical_pattern_recon(prtrks[ht])
    tracker = hl_prtr["tracker"]
    detector = "upstream" if tracker == 0 else "downstream"
    container[detector].append(hl_prtr)
  return container

#########################################################################################
  # Collects information from track points and tracks.  
def Tracks(tracks):
  container = {"upstream":[], "downstream":[]}
  for tk in range(len(tracks)):
    track = fill_tracks(tracks[tk])
    tracker = track["tracker"]
    detector = "upstream" if tracker == 0 else "downstream"
    container[detector].append(track)
  return container

#########################################################################################
  # Collects information from TOF space points.
def TOF(tofpts):
  container = []
  for tf in range(len(tofpts)):
    tofpt = fill_tof(tofpts[tf])
    container.append(tofpt)
  return container

#########################################################################################
  # Collects information from MC Virtual Hits
def Virtual_Hits(virt):
  container = {"upstream":{1:[], 2:[], 3:[], 4:[], 5:[]}, \
             "downstream":{1:[], 2:[], 3:[], 4:[], 5:[]}, \
             "uncat":[]}
  for vi in range(len(virt)):
    virtual  = fill_virtual_hits(virt[vi])
    detector = virtual["detector"]
    if detector == "uncat":
      container[detector].append(virtual)
    else:
      station  = virtual["station"]
      container[detector][station].append(virtual)
  return container


def fill_tof(tofpt):
  temp               = {}
  temp["TOF"]        = tofpt.GetDetector()
  temp["time"]       = tofpt.GetTime()
  temp["x_pos"]      = tofpt.GetGlobalPosX()
  temp["y_pos"]      = tofpt.GetGlobalPosY()
  temp["z_pos"]      = tofpt.GetGlobalPosZ()
  temp["x_err"]      = tofpt.GetGlobalPosXErr()
  temp["y_err"]      = tofpt.GetGlobalPosYErr()
  temp["z_err"]      = tofpt.GetGlobalPosZErr()
  temp["dt"]         = tofpt.GetDt()
  temp["part_event"] = tofpt.GetPartEventNumber()
  temp["phys_event"] = tofpt.GetPhysEventNumber()
  temp["x_slab"]     = tofpt.GetSlabx()
  temp["y_slab"]     = tofpt.GetSlaby()
  return temp

def fill_straight_pattern_recon(prtrk):
  temp               = {}
  temp["tracker"]    = prtrk.get_tracker()
  temp["z_mom"]      = prtrk.get_reference_momentum().z()
  temp["x_mom"]      = prtrk.get_reference_momentum().x()
  temp["y_mom"]      = prtrk.get_reference_momentum().y()
  temp["x_0"]        = prtrk.get_x0()
  temp["y_0"]        = prtrk.get_y0()
  temp["m_x"]        = prtrk.get_mx()
  temp["m_y"]        = prtrk.get_my()
  temp["x_chisq"]    = prtrk.get_x_chisq()
  temp["y_chisq"]    = prtrk.get_y_chisq()
  temp["num_pnt"]    = prtrk.get_num_points()
  temp["triplets"]   = prtrk.get_num_triplets()
  temp["chi_sq"]     = prtrk.get_chi_squared()
  temp["ref_pos_x"]  = prtrk.get_reference_position().x()
  temp["ref_pos_y"]  = prtrk.get_reference_position().y()
  temp["ref_pos_z"]  = prtrk.get_reference_position().z()
  temp["ndf"]        = prtrk.get_ndf()
  temp["type"]       = prtrk.get_type()
  temp["covariance"] = prtrk.get_covariance()
  temp["seeds"]      = []
  temp["pr_covariance"] = prtrk.get_covariance()
  seed = Space_Points(prtrk.get_spacepoints())
  temp["seeds"].append(seed)
  return temp

def fill_space_points(space):
  temp               = {}
  temp["tracker"]    = space.get_tracker()
  temp["pe"]         = space.get_npe()
  temp["station"]    = space.get_station()
  temp["spill"]      = space.get_spill()
  temp["event"]      = space.get_event()
  temp["x_pos"]      = space.get_position().x()
  temp["y_pos"]      = space.get_position().y()
  temp["z_pos"]      = space.get_position().z()
  temp["time"]       = space.get_time()
  temp["chi_sq"]     = space.get_chi2()
  temp["type"]       = 3 if space.get_type() == "triplet" else 2
  temp["x_glob_pos"] = space.get_global_position().x()
  temp["y_glob_pos"] = space.get_global_position().y()
  temp["z_glob_pos"] = space.get_global_position().z()
  temp["time_err"]   = space.get_time_error()
  temp["time_res"]   = space.get_time_res()
  temp["clusters"]   = []
  cluster = Clusters(space.get_channels())
  temp["clusters"].append(cluster)
  return temp

def fill_clusters(cluster):
  temp                = {}
  temp["spill"]       = cluster.get_spill()
  temp["event"]       = cluster.get_event()
  temp["tracker"]     = cluster.get_tracker()
  temp["station"]     = cluster.get_station()
  temp["plane"]       = cluster.get_plane()
  temp["pe"]          = cluster.get_npe()
  temp["time"]        = cluster.get_time()
  temp["channel"]     = cluster.get_channel()
  temp["direction_x"] = cluster.get_direction().x()
  temp["direction_y"] = cluster.get_direction().y()
  temp["direction_z"] = cluster.get_direction().z()
  temp["pos_x"]       = cluster.get_position().x()
  temp["pos_y"]       = cluster.get_position().y()
  temp["pos_z"]       = cluster.get_position().z()
  temp["alpha"]       = cluster.get_alpha()
  temp["type"]        = len(cluster.get_digits())
  temp["digits"]      = []
  digit = Digits(cluster.get_digits())
  temp["digits"].append(digit)
  return temp

def fill_digits(digit):
  temp            = {}
  temp["spill"]   = digit.get_spill()
  temp["event"]   = digit.get_event()
  temp["tracker"] = digit.get_tracker()
  temp["station"] = digit.get_station()
  temp["plane"]   = digit.get_plane()
  temp["pe"]      = digit.get_npe()
  temp["channel"] = digit.get_channel()
  temp["time"]    = digit.get_time()
  temp["adc"]     = digit.get_adc()
  return temp

def fill_helical_pattern_recon(prtrk):
  temp                  = {}
  temp["phi"]           = prtrk.get_phi()
  temp["tracker"]       = prtrk.get_tracker()
  temp["charge"]        = prtrk.get_charge()
  temp["pos0_x"]        = prtrk.get_pos0().x()
  temp["pos0_y"]        = prtrk.get_pos0().y()
  temp["pos0_z"]        = prtrk.get_pos0().z()
  temp["phi0"]          = prtrk.get_phi0()
  temp["dsdz"]          = prtrk.get_dsdz()
  temp["radius"]        = prtrk.get_R()
  temp["line_sz_c"]     = prtrk.get_line_sz_c()
  temp["line_sz_chisq"] = prtrk.get_line_sz_chisq()
  temp["circle_x0"]     = prtrk.get_circle_x0()
  temp["circle_y0"]     = prtrk.get_circle_y0()
  temp["circle_chisq"]  = prtrk.get_circle_chisq()
  temp["point_spread"]  = prtrk.get_point_spread()
  temp["type"]          = prtrk.get_type()
  temp["covariance"]    = prtrk.get_covariance()
  temp["ndf"]           = prtrk.get_ndf()
  temp["num_pnt"]       = prtrk.get_num_points()
  temp["triplets"]      = prtrk.get_num_triplets()
  temp["ref_pos_x"]     = prtrk.get_reference_position().x()
  temp["ref_pos_y"]     = prtrk.get_reference_position().y()
  temp["ref_pos_z"]     = prtrk.get_reference_position().z()
  temp["z_mom"]         = prtrk.get_reference_momentum().z()
  temp["x_mom"]         = prtrk.get_reference_momentum().x()
  temp["y_mom"]         = prtrk.get_reference_momentum().y()
  temp["chi_sq"]        = prtrk.get_chi_squared()
  temp["seeds"]         = []
  seed = Space_Points(prtrk.get_spacepoints())
  temp["seeds"].append(seed)
  return temp

def fill_tracks(track):
  temp                 = {}
  num_data             = 0
  temp["tracker"]      = track.tracker()
  temp["chi_sq"]       = track.chi2()
  temp["ndf"]          = track.ndf()
  temp["p_value"]      = track.P_value()
  temp["charge"]       = track.charge()
  temp["algorithm"]    = track.GetAlgorithmUsed()
  temp["rating"]       = track.GetRating()
  temp["good"]         = track.IsGood()
  temp["prtrks"]       = []
  temp["track_points"] = []
  if temp["algorithm"] == 0:
    prtrk = fill_straight_pattern_recon(track.pr_track_pointer_straight())
    temp["prtrks"].append(prtrk)
  else:
    prtrk = fill_helical_pattern_recon(track.pr_track_pointer_helical()[pr])
    temp["prtrks"].append(prtrk)
  for tp in range(len(track.scifitrackpoints())):
    trpint                     = track.scifitrackpoints()[tp]
    temp2                      = {}
    temp2["tracker"]           = trpint.tracker()
    temp2["station"]           = trpint.station()
    temp2["plane"]             = trpint.plane()
    temp2["channel"]           = trpint.channel()
    temp2["chi_sq"]            = trpint.chi2()
    temp2["x_pos"]             = trpint.pos().x()
    temp2["y_pos"]             = trpint.pos().y()
    temp2["z_pos"]             = trpint.pos().z()
    temp2["x_mom"]             = trpint.mom().x()
    temp2["y_mom"]             = trpint.mom().y()
    temp2["z_mom"]             = trpint.mom().z()
    temp2["covariance"]        = trpint.covariance()
    temp2["errors"]            = trpint.errors()
    temp2["pull"]              = trpint.pull()
    temp2["residual"]          = trpint.residual()
    temp2["smoothed_residual"] = trpint.smoothed_residual()
    temp2["spill"]             = trpint.spill()
    temp2["event"]             = trpint.event()
    temp2["has_data"]          = trpint.has_data()
    if temp2["has_data"] == True:
      num_data =+ 1
    temp["track_points"].append(temp2)
  temp["number_data_points"] = num_data
  return temp

def fill_virtual_hits(virtual):
  temp = {}
  temp["station_id"]  = virtual.GetStationId()
  temp["track_id"]    = virtual.GetTrackId()
  temp["time"]        = virtual.GetTime()
  temp["mass"]        = virtual.GetMass()
  temp["charge"]      = virtual.GetCharge()
  temp["proper_time"] = virtual.GetProperTime()
  temp["path_length"] = virtual.GetPathLength()
  temp["spin"]        = virtual.GetSpin
  temp["x_pos"]       = virtual.GetPosition().x()
  temp["y_pos"]       = virtual.GetPosition().y()
  temp["z_pos"]       = virtual.GetPosition().z()
  temp["x_B"]         = virtual.GetBField().x()
  temp["y_B"]         = virtual.GetBField().y()
  temp["z_B"]         = virtual.GetBField().z()
  temp["x_E"]         = virtual.GetBField().x()
  temp["y_E"]         = virtual.GetBField().y()
  temp["z_E"]         = virtual.GetBField().z()
  temp["x_mom"]       = virtual.GetMomentum().x()
  temp["y_mom"]       = virtual.GetMomentum().y()
  temp["z_mom"]       = virtual.GetMomentum().z()
  temp["part_id"]     = virtual.GetParticleId()
  temp["detector"]    = "uncat"
  virt_map = _config.fill_config["Virtual_to_Tracker_Map"]
  if temp["station_id"] in virt_map:
    temp["tracker"]  = virt_map[temp["station_id"]][0]
    temp["station"]  = virt_map[temp["station_id"]][1]
    temp["plane"]    = virt_map[temp["station_id"]][2]
    temp["detector"] = "upstream" if temp["tracker"] == 0 else "downstream"
  return temp