#!/usr/bin/env python
import libMausCpp
import ROOT
from ROOT import gROOT
import os
from os import listdir
import _ctypes
from array import *
import argparse
from config_maus_analysis import check_config as _config
import fill_maus_analysis as _fill
import analyize_maus_analysis as _analyize

class Process:
  def __init__ (self):
    self.a = _analyize.Analyize()
#    self.a = _analyize.Analyize()
#    self.a.Initialize()
    self.Read_Spills()
    # self.Output()

#########################################################################################
  # Reads MICE spill data
  def Read_Spills(self):
  # Creates the container ROOT files and loads up the MAUS processed file.
    self.Make_ROOT()
    file_in = self.Load_file()
    print "Reading MAUS processed file: ",file_in
    root_file = ROOT.TFile(file_in, "READ") 
  # Checks spill/event is good data
    data = ROOT.MAUS.Data()
    tree = root_file.Get("Spill")
    if not tree:
      return
    tree.SetBranchAddress("data", data)
    peat_count = 0
    event_cut = _config["event_cut"]
    if event_cut > tree.GetEntries() or event_cut <= 0:
      event_cut = tree.GetEntries()
    for i in range(event_cut):
      peat_count += 1
      if (peat_count % _config["event_out"] == 0):
        print "Filling event: ",peat_count, "/", event_cut
      tree.GetEntry(i)
      self.spill = data.GetSpill()
      if not self.spill.GetDaqEventType() == "physics_event":
        continue
      if self.spill.GetReconEvents().size() == 0:
        print "No Recon Events"
        continue

  # Fills the ROOT containers with data from the MAUS files
      self.Process_Event()
    self.Output()

#########################################################################################
  # Reads a single MICE event and sets up all the contianers that will be needed
  #   to store the data for that event. Calls functions to read each type of 
  #   MICE data to fill those contianers. Finally calls a function to start 
  #   analysis of the event using collected data and flags.  
  def Process_Event(self):
    for rc in range(self.spill.GetReconEvents().size()):
      _config["counter"]["total_events"]+=1
      self.is_recon = False
      self.is_mc    = False
      if not self.spill.GetReconEvents()[rc] == 0:
        self.is_recon = True
        recon_event   = self.spill.GetReconEvents()[rc]
      if not self.spill.GetMCEvents().size() == 0:
        self.is_mc    = True
        mc_event      = self.spill.GetMCEvents()[rc]
      self.data = _fill.Fill_from_Data(recon_event, mc_event, \
                                       self.is_recon, self.is_mc)
      self.Call_Analysis()

#########################################################################################
  # Start Analysis Functions
#########################################################################################

#########################################################################################
  # Generates a new copy of the virtual hits to tracker station/plane map
  #  from comparison against track points.  It is probably a good idea to run
  #  over at least 10 spills to make sure each station is activated.
  # The output file can be copied and pasted in config_math.py in the function
  #   fill_config and into the Virtual_to_Tracker_Map variable
  def Generate_Virtual_Map(self, virtuals, tracks):
    config_dict = {}
    for detector in virtuals:
      for event in virtuals[detector]:
        for virtual in virtuals["uncat"][event]:
          for detector in tracks:
            for track in tracks[detector][event]:
              for space in track["track_points"]:
                if abs(virtual["z_pos"] - space["z_pos"]) < 0.012:
                  config_dict[virtual["station_id"]] = [space["tracker"], \
                                                        space["station"], \
                                                        space["plane"]]
    file = open('virtual_map', 'w')
    file.write(config_dict)

#########################################################################################
  #
  def SP_to_Virt(self, virtuals, spaces):
    for detector in spaces:
      for station in spaces[detector]:
        if len(spaces[detector][station]) == 1:
          for virtual in virtuals[detector][station]:
            if virtual["plane"] == 0:
              space   = spaces[detector][station][0]
              x_res   = virtual["x_pos"] - space["x_glob_pos"]
              y_res   = virtual["y_pos"] - space["y_glob_pos"]
              i = "U" if detector == "upstream" else "D"
              self.SP_to_Virt_Hist[i][station].Fill(x_res, y_res)
              self.SP_to_Virt_Hist_x[i][station].Fill(x_res)
              self.SP_to_Virt_Hist_y[i][station].Fill(y_res)
            
#########################################################################################
  #
  def SP_Fill_ROOT(self, space_points):
    for detector in space_points:
      for station in space_points[detector]:
        if len(space_points[detector][station]) == 1:
          x_pos   = space_points[detector][station][0]["x_pos"]
          y_pos   = space_points[detector][station][0]["y_pos"]
          i = "U" if detector == "upstream" else "D"
          self.SP_Pos[i][station].Fill(x_pos, y_pos)
          if space_points[detector][station][0]["type"] == 3:
            self.Triplet_Pos[i][station].Fill(x_pos, y_pos)
          else:
            self.Doublet_Pos[i][station].Fill(x_pos, y_pos)

#########################################################################################
  #
  def Virt_Fill_ROOT(self, virtuals):
    for detector in virtuals:
      if detector == "uncat":
        for virtual in virtuals[detector]:
          x_pos = virtual["x_pos"]
          y_pos = virtual["y_pos"]
          self.Virt_Pos["UnCat"].Fill(x_pos, y_pos)
      else:
        for station in virtuals[detector]:
          for virtual in virtuals[detector][station]:
            x_pos   = virtual["x_pos"]
            y_pos   = virtual["y_pos"]
            i = "U" if detector == "upstream" else "D"
            self.Virt_Pos[i][station].Fill(x_pos, y_pos)

#########################################################################################
  #
  def TOF_Timing_Info(self, up_tof, down_tof):
    if len(up_tof) == 1 and \
       len(down_tof) == 1:
      timing = up_tof["time"] - down_tof["time"]
      return timing
    else:
      timing = -10000.
      return timing
      
#########################################################################################
  # Reads in TOF1 and TOF2 positions and draws a line between the two.  Space
  #   points are transformed into global coordinates and the residuals between
  #   where the TOF to TOF line cross the tracker plane and space points are 
  #   calculated.
  def TOF_to_TOF_Tkr_Res(self, tracks, tof1, tof2):
    for detector in tracks
      if tracks[detector]["triples"] < _config["TtT_trip_req"]:
        continue
      TOF_distance = tof1["z_pos"] - tof2["z_pos"]
      x_change = tof1["x_pos"] - tof2["x_pos"]
      y_change = tof1["y_pos"] - tof2["y_pos"]
      x_slope  = x_change/TOF_distance
      y_slope  = y_change/TOF_distance

      for spaces in track[detector]["seeds"]:
        for station in spaces[detector]:
          space = spaces[detector][station]
          z_pos = space["z_glob_pos"]
          distance = z_pos - tof1["z_pos"]
          expected_x = tof1["x_pos"] + distance * x_slope
          residual_x = space["x_glob_pos"] - expected_x
          expected_y = tof1["y_pos"] + distance * y_slope
          residual_y = space["y_glob_pos"] - expected_y
          self.tof_to_tof_residual["up"][st].Fill(residual_x,residual_y)
          self.tof_to_tof_residual_x["up"][st].Fill(residual_x)
          self.tof_to_tof_residual_y["up"][st].Fill(residual_y)

#########################################################################################
  # Backend functions
#########################################################################################

#########################################################################################
  # Calls cuts and analysis routines.  Passes in data containers, event numbers
  #   and histograms to be filled.
  def Call_Analysis(self):
    TOF_timing = {"upstream":False, "downstream":False}
    if _config["ignore_Generate_Virtual_Map"] == False:
      self.Generate_Virtual_Map(self.data["virtual_points"], \
                                self.data["tracker_tracks"])

    if _config["ignore_SP_to_Virt"] == False and \
               "virtual_points" in self.data and \
               "tracker_space_points" in self.data:
      self.SP_to_Virt(self.data["virtual_points"], \
                           self.data["tracker_space_points"])

    if _config["ignore_SP_Fill_ROOT"] == False and \
               "tracker_space_points" in self.data:
      self.SP_Fill_ROOT(self.data["tracker_space_points"])

    if _config["ignore_Virt_Fill_ROOT"] == False and \
               "virtual_points" in self.data:
      self.Virt_Fill_ROOT(self.data["virtual_points"])

    if _config["ignore_TOF_Timing_Info"] == False and \
               "TOF0_space_points" in self.data and \
               "TOF1_space_points" in self.data:
      time = self.TOF_Timing_Info(self.data["TOF0_space_points"], \
                                  self.data["TOF1_space_points"])
      if time > _config["upstream_Tmin"] and \
         time < _config["upstream_Tmax"]:
        TOF_timing["upstream"] = True
    if _config["ignore_TOF_Timing_Info"] == False and \
               "TOF1_space_points" in self.data and \
               "TOF2_space_points" in self.data:
      time = self.TOF_Timing_Info(self.data["TOF1_space_points"], \
                                  self.data["TOF2_space_points"])
      if time > _config["Downstream_Tmin"] and \
         time < _config["Downstream_Tmax"]:
        TOF_timing["downstream"] = True

    if _config["ignore_Station_Alignment"] == False and \
               "tracker_straight_pr" in self.data:
      self.a.StS_Collect_Space_Points(self.data["tracker_straight_pr"])
      
    if _config["ignore_TOF_to_TOF_Tkr_Res"]  == False and \
               "tracker_straight_pr" in self.data and \
               "TOF1_space_points" in self.data and \
               "TOF2_space_points" in self.data and \
               TOF_timing["downstream"] == True:
      self.TOF_to_TOF_Tkr_Res(self.data["tracker_straight_pr"] \
                              self.data["TOF1_space_points"] \
                              self.data["TOF2_space_points"])

#########################################################################################
  #
  def Output(self):
    self.a.Station_Alignment()
    out_root = ROOT.TFile(_config["output_file"],'RECREATE')
    self.Virt_Pos["UnCat"].Write()
    for st in range(1, 6):
      self.SP_to_Virt_Hist["U"][st].Write()
      self.SP_to_Virt_Hist["D"][st].Write()
      self.SP_to_Virt_Hist_x["U"][st].Write()
      self.SP_to_Virt_Hist_x["D"][st].Write()
      self.SP_to_Virt_Hist_y["U"][st].Write()
      self.SP_to_Virt_Hist_y["D"][st].Write()
      self.Triplet_Pos["U"][st].Write()
      self.Triplet_Pos["D"][st].Write()
      self.Doublet_Pos["U"][st].Write()
      self.Doublet_Pos["D"][st].Write()
      self.SP_Pos["U"][st].Write()
      self.SP_Pos["D"][st].Write()
#      self.Virt_Pos["U"][st].Write()
#      self.Virt_Pos["D"][st].Write()
      
    gROOT.SetBatch()
    c1 = ROOT.TCanvas("c1","Space Point Positions UpStream",1200, 800)
    ROOT.SetOwnership(c1, False)
    c1.Divide(3,2)
    c1.cd(1)
    self.SP_to_Virt_Hist["U"][1].Draw()
    c1.cd(2)
    self.SP_to_Virt_Hist["U"][2].Draw()
    c1.cd(3)
    self.SP_to_Virt_Hist["U"][3].Draw()
    c1.cd(4)
    self.SP_to_Virt_Hist["U"][4].Draw()
    c1.cd(5)
    self.SP_to_Virt_Hist["U"][5].Draw()
    c1.Print("Recon_Truth_Residual_Up.pdf")

    c2 = ROOT.TCanvas("c2","Space Point Positions DownStream",1200, 800)
    ROOT.SetOwnership(c2, False)
    c2.Divide(3,2)
    c2.cd(1)
    self.SP_to_Virt_Hist["D"][1].Draw()
    c2.cd(2)
    self.SP_to_Virt_Hist["D"][2].Draw()
    c2.cd(3)
    self.SP_to_Virt_Hist["D"][3].Draw()
    c2.cd(4)
    self.SP_to_Virt_Hist["D"][4].Draw()
    c2.cd(5)
    self.SP_to_Virt_Hist["D"][5].Draw()
    c2.Print("Recon_Truth_Residual_Down.pdf")
    
    c3 = ROOT.TCanvas("c3","Virtual Point Positions UpStream",1200, 800)
    ROOT.SetOwnership(c3, False)
    c3.Divide(3,2)
    c3.cd(1)
    self.Virt_Pos["U"][1].Draw()
    c3.cd(2)
    self.Virt_Pos["U"][2].Draw()
    c3.cd(3)
    self.Virt_Pos["U"][3].Draw()
    c3.cd(4)
    self.Virt_Pos["U"][4].Draw()
    c3.cd(5)
    self.Virt_Pos["U"][5].Draw()

    c4 = ROOT.TCanvas("c4","Virtual Point Positions DownStream",1200, 800)
    ROOT.SetOwnership(c4, False)
    c4.Divide(3,2)
    c4.cd(1)
    self.Virt_Pos["D"][1].Draw()
    c4.cd(2)
    self.Virt_Pos["D"][2].Draw()
    c4.cd(3)
    self.Virt_Pos["D"][3].Draw()
    c4.cd(4)
    self.Virt_Pos["D"][4].Draw()
    c4.cd(5)
    self.Virt_Pos["D"][5].Draw()

    out_root.Close()
    raw_input("Press Enter to Exit")

#########################################################################################
  # Creates the container root file that will be passed around script.
  def Make_ROOT(self):
    print "Creating empty ROOT file"
    self.SP_to_Virt_Hist   = {"U":{},"D":{}}
    self.SP_to_Virt_Hist_x = {"U":{},"D":{}}
    self.SP_to_Virt_Hist_y = {"U":{},"D":{}}
    self.SP_Pos            = {"U":{},"D":{}}
    self.Virt_Pos          = {"U":{},"D":{}}
    self.Triplet_Pos       = {"U":{},"D":{}}
    self.Doublet_Pos       = {"U":{},"D":{}}
    self.Virt_Pos["UnCat"] = ROOT.TH2D("Virt_Pos", \
                                         "Virtual Space Point Position", \
                                         400, -400 , 400, 400, -400 , 400)
    for st in range(1,6):
      self.SP_to_Virt_Hist["U"][st] = ROOT.TH2D("SP_to_Virt_U%i" % st, \
                                                "Residual Virtual Plane to \
                                                Tracker Space Point TKU %i " % st, \
                                                500, -10, 10, 500, -10, 10)
      self.SP_to_Virt_Hist["D"][st] = ROOT.TH2D("SP_to_Virt_D%i" % st, \
                                                "Residual Virtual Plane to \
                                                Tracker Space Point TKD %i " % st, \
                                                500, -10, 10, 500, -10, 10)
      self.SP_to_Virt_Hist_x["U"][st] = ROOT.TH1D("SP_to_Virt_U%i_X" % st, \
                                                  "Residual Virtual Plane to \
                                                  Tracker Space Point TKU %i X" % st, \
                                                  500, -10, 10)
      self.SP_to_Virt_Hist_x["D"][st] = ROOT.TH1D("SP_to_Virt_D%i_X" % st, \
                                                  "Residual Virtual Plane to \
                                                  Tracker Space Point TKD %i X" % st, \
                                                  500, -10, 10)
      self.SP_to_Virt_Hist_y["U"][st] = ROOT.TH1D("SP_to_Virt_U%i_Y" % st, \
                                                  "Residual Virtual Plane to \
                                                  Tracker Space Point TKU %i Y" % st, \
                                                  500, -10, 10)
      self.SP_to_Virt_Hist_y["D"][st] = ROOT.TH1D("SP_to_Virt_D%i_Y" % st, \
                                                  "Residual Virtual Plane to \
                                                  Tracker Space Point TKD %i Y" % st, \
                                                  500, -10, 10)
      self.Triplet_Pos["U"][st]   = ROOT.TH2D("Triplet_Pos_U%i" %st, \
                                         "Tracker SP Triplets TKU %i" %st, \
                                         250, -250 , 250, 250, -250 , 250)
      self.Triplet_Pos["D"][st]   = ROOT.TH2D("Triplet_Pos_D%i" %st, \
                                         "Tracker SP Triplets TKD %i" %st, \
                                         250, -250 , 250, 250, -250 , 250)
      self.Doublet_Pos["U"][st]   = ROOT.TH2D("Doublet_Pos_U%i" %st, \
                                         "Tracker SP Doublets TKU %i" %st, \
                                         250, -250 , 250, 250, -250 , 250)
      self.Doublet_Pos["D"][st]   = ROOT.TH2D("Doublet_Pos_D%i" %st, \
                                         "Tracker SP Doublets TKD %i" %st, \
                                         250, -250 , 250, 250, -250 , 250)
      self.SP_Pos["U"][st]   = ROOT.TH2D("SP_Pos_U%i" %st, \
                                         "Tracker Space Point Position TKU %i" %st, \
                                         250, -250 , 250, 250, -250 , 250)
      self.SP_Pos["D"][st]   = ROOT.TH2D("SP_Pos_D%i" %st, \
                                         "Tracker Space Point Position TKD %i" %st, \
                                         250, -250 , 250, 250, -250 , 250)
      self.Virt_Pos["U"][st] = ROOT.TH2D("Virt_Pos_U%i" %st, \
                                         "Virtual Space Point Position TKU %i" %st, \
                                         250, -250 , 250, 250, -250 , 250)
      self.Virt_Pos["D"][st] = ROOT.TH2D("Virt_Pos_D%i" %st, \
                                         "Virtual Space Point Position TKD %i" %st, \
                                         250, -250 , 250, 250, -250 , 250)

#########################################################################################
  # Searches predefined data directory to find specified processed MAUS file.
  def Load_file(self):
    for input_file in listdir(_config["data_directory"]):
      if _config["data_identifier"] in input_file and ".root" in input_file:
        file_in = _config["data_directory"] + input_file
    return file_in

#########################################################################################
if __name__ == "__main__":
  Process()