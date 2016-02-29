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
import analyize_maus_analysis as _analysis
import station_alignment as _st_align

class Process:
  def __init__ (self):
    self.st_align = _st_align.ST_Alignment()
    self.analysis = _analysis.Analysis()
    self.Read_Spills()

#########################################################################################
  # Reads MICE spill data
  def Read_Spills(self):
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
  # Reads in TOF1 and TOF2 positions and draws a line between the two.  Space
  #   points are transformed into global coordinates and the residuals between
  #   where the TOF to TOF line cross the tracker plane and space points are 
  #   calculated.
  def TOF_to_TOF_Tkr_Res(self, tracks, tof1, tof2):
    for detector in tracks:
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
      self.analysis.SP_to_Virt(self.data["virtual_points"], \
                               self.data["tracker_space_points"])

    if _config["ignore_SP_Fill_ROOT"] == False and \
               "tracker_space_points" in self.data:
      self.analysis.SP_Fill_ROOT(self.data["tracker_space_points"])

    if _config["ignore_Virt_Fill_ROOT"] == False and \
               "virtual_points" in self.data:
      self.analysis.Virt_Fill_ROOT(self.data["virtual_points"])

    if _config["ignore_TOF_Timing_Info"] == False and \
               "TOF0_space_points" in self.data and \
               "TOF1_space_points" in self.data:
      time = self.analysis.TOF_Timing_Info(self.data["TOF0_space_points"], \
                                           self.data["TOF1_space_points"])
      if time > _config["upstream_Tmin"] and \
         time < _config["upstream_Tmax"]:
        TOF_timing["upstream"] = True
    if _config["ignore_TOF_Timing_Info"] == False and \
               "TOF1_space_points" in self.data and \
               "TOF2_space_points" in self.data:
      time = self.analysis.TOF_Timing_Info(self.data["TOF1_space_points"], \
                                           self.data["TOF2_space_points"])
      if time > _config["Downstream_Tmin"] and \
         time < _config["Downstream_Tmax"]:
        TOF_timing["downstream"] = True

    if _config["ignore_Station_Alignment"] == False and \
               "tracker_straight_pr" in self.data:
      self.st_align.StS_Collect_Space_Points(self.data["tracker_straight_pr"])
      
    if _config["ignore_TOF_to_TOF_Tkr_Res"]  == False and \
               "tracker_straight_pr" in self.data and \
               "TOF1_space_points" in self.data and \
               "TOF2_space_points" in self.data and \
               TOF_timing["downstream"] == True:
      self.TOF_to_TOF_Tkr_Res(self.data["tracker_straight_pr"], \
                              self.data["TOF1_space_points"],   \
                              self.data["TOF2_space_points"])

#########################################################################################
  #
  def Output(self):
    self.st_align.Station_Alignment()
    self.analysis.Write()
    raw_input("Press Enter to Exit")

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