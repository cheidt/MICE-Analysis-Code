#!/usr/bin/env python
import libMausCpp
import ROOT
from ROOT import gROOT

#import os
from os import listdir
import _ctypes
from array import *
import argparse
import subprocess

from config_maus_analysis import check_config as _config
import fill_maus_analysis as _fill
import analyize_maus_analysis as _analysis
import station_alignment as _st_align
import output_maus_analysis as _output

class Process:
  def __init__ (self):
    self.st_align = _st_align.ST_Alignment()
    self.analysis = _analysis.Analysis()
    self.Read_Spills()

#########################################################################################
  # Reads MICE spill data
  def Read_Spills(self):
    _output.Message("INITALIZING LOG FILE: ",_config["log_file"])
    clear_log = open(_config["log_file"],"w")
    clear_log.close()
    file_in = self.Load_file()
    _output.Message("READING MAUS FILE: ", file_in, exc=True)
    root_file = ROOT.TFile(file_in, "READ") 
  # Checks spill/event is good data
    data = ROOT.MAUS.Data()
    tree = root_file.Get("Spill")
    if not tree:
      return
    tree.SetBranchAddress("data", data)
    event_cut = _config["event_cut"]
    min_event = _config["min_event"]
    if event_cut > tree.GetEntries() or event_cut <= 0:
      event_cut = tree.GetEntries()
    for i in range(min_event, event_cut):
      if (i % _config["event_out"] == 0):
        _output.Message("Filling event: ", i, "/", event_cut, exc=True)
      tree.GetEntry(i)
      self.spill = data.GetSpill()
      if not self.spill.GetDaqEventType() == "physics_event":
        continue
      if self.spill.GetReconEvents().size() == 0 and \
         self.spill.GetMCEvents().size() == 0:
        _output.Message("No data in event")
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
    if not self.spill.GetMCEvents().size() == 0:
      for rc in range(self.spill.GetMCEvents().size()):
        _config["counter"]["total_events"]+=1
        if self.spill.GetReconEvents()[rc]:
          self.data = _fill.Fill_from_Data( \
                            recon = self.spill.GetReconEvents()[rc], \
                            mc    = self.spill.GetMCEvents()[rc])
          self.Call_Analysis()
        else:
          self.data = _fill.Fill_from_Data( \
                            mc    = self.spill.GetMCEvents()[rc])
          self.Call_Analysis()
    else:
      for rc in range(self.spill.GetReconEvents().size()):
        _config["counter"]["total_events"]+=1
        self.data = _fill.Fill_from_Data( \
                          recon = self.spill.GetReconEvents()[rc])
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

    time = -10000
    if _config["ignore_TOF_Timing_Info"] == False and \
               "TOF0_space_points" in self.data and \
               "TOF1_space_points" in self.data:
      time = self.analysis.TOF_Timing_Info(self.data["TOF0_space_points"], \
                                           self.data["TOF1_space_points"])
      if time > _config["upstream_Tmin"] and \
         time < _config["upstream_Tmax"]:
        TOF_timing["upstream"] = True
    time = -10000
    if _config["ignore_TOF_Timing_Info"] == False and \
               "TOF1_space_points" in self.data and \
               "TOF2_space_points" in self.data:
      time = self.analysis.TOF_Timing_Info(self.data["TOF1_space_points"], \
                                           self.data["TOF2_space_points"])
      if time > _config["downstream_Tmin"] and \
         time < _config["downstream_Tmax"]:
        TOF_timing["downstream"] = True

    if _config["ignore_Station_Alignment"] == False and \
               "tracker_straight_pr" in self.data:
      self.st_align.StS_Collect_Space_Points(self.data["tracker_straight_pr"])

    if _config["ignore_TOF_Tkr_Res"]  == False and \
               "tracker_tracks" in self.data and \
               "TOF1_space_points" in self.data and \
               "TOF2_space_points" in self.data and \
               TOF_timing["downstream"] == True:
      self.analysis.TOF_Tk_Res(self.data["tracker_tracks"], \
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
      if _config["file_name"] in input_file and ".root" in input_file:
        file_in = _config["data_directory"] + input_file
    return file_in

#########################################################################################
if __name__ == "__main__":
  Process()