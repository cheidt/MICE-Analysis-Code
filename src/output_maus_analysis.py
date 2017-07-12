import ROOT
import os
import json
import pprint
from config_maus_analysis import out_config as _config

class Output(object):
  def __init__(self, type):
    self.name = type
    self.hist_container = {}

  def Fill(self, name, title, *args, **kwargs):
    if "detector" in kwargs:
      title = title + " " + kwargs["detector"]
      name  = name  + "_" + kwargs["detector"]
    if "station" in kwargs:
      title = title + " " + str(kwargs["station"])
      name  = name  + "_" + str(kwargs["station"])
    if "plane" in kwargs:
      title = title + " " + str(kwargs["plane"])
      name  = name  + "_" + str(kwargs["plane"])
    if "iteration" in kwargs:
      title = title + " itr " + str(kwargs["iteration"])
      name  = name  + "_"     + str(kwargs["iteration"])
    if "channel" in kwargs:
      title = title + " " + str(kwargs["channel"])
      name  = name  + "_" + str(kwargs["channel"])
    if "weight" in kwargs:
      w = kwargs["weight"]
    else:
      w = 1.0
    if name not in self.hist_container:
      self.Initialize(name, title, args)
    else:
      if len(args) == 8:
        self.hist_container[name].Fill(args[0], args[1], w)
      else:
        self.hist_container[name].Fill(args[0], w)

  def Graph(self, name, title, *args, **kwargs):
    if name not in self.hist_container:
      self.hist_container[name] = ROOT.TGraph()
      self.hist_container[name].SetName(name)
      self.hist_container[name].SetTitle(title)
      self.hist_container[name].SetMarkerStyle(3)
      self.hist_container[name].SetLineColor(0)
      self.hist_container[name].SetLineColorAlpha(0, 1.0)
      if "color" in kwargs:
        self.hist_container[name].SetMarkerColor(kwargs["color"])
    self.hist_container[name].SetPoint(self.hist_container[name].GetN(), args[0], args[1])

  def Initialize(self, name, title, args):
    if len(args) == 8:
      self.hist_container[name] = ROOT.TH2D(name, title, args[2], args[3], \
                                                         args[4], args[5], \
                                                         args[6], args[7])
      self.hist_container[name].SetDirectory(0)
    else:
      self.hist_container[name] = ROOT.TH1D(name, title,args[1], args[2], \
                                                        args[3])
      self.hist_container[name].SetDirectory(0)

  def Add_Graphs(self, name, title, *args):
    if name not in self.hist_container:
      self.hist_container[name] = ROOT.TMultiGraph(name, title)
    for graph in args:
      self.hist_container[name].Add(self.hist_container[graph])

  def Write(self):
    file_name = _config["output_dir"] + self.name + _config["output_file"]
    print "WRITTING ANALYSIS FILE: ",file_name
    out_root = ROOT.TFile(file_name,'RECREATE')
    for name in self.hist_container:
      self.hist_container[name].Write()
    out_root.Close()

  def Draw(self, name, *args):
    self.hist_container[name].Clear()
    for i in range(len(args)):
      self.hist_container[name].SetOption(args[i])
    self.hist_container[name].Draw()

  def Gaus_Fit(self, low, high, name):
    g1 = ROOT.TF1('m1','gaus', low, high)
    self.hist_container[name].Fit(g1,'R')

  def Delete(self, name):
    if name in self.hist_container:
      del self.hist_container[name]
    else:
      pass

  def Fit(self):
    for name in self.hist_container:
      g4 = ROOT.TF1('m1','landau',22,27)
      g5 = ROOT.TF1('m2','landau',27,34)
      g6 = ROOT.TF1('m3','landau',32,40)
      #total = ROOT.TF1('mstotal','landau(0)+landau(3)+landau(6)',20,50)

      g1 = ROOT.TF1('m1','gaus',22,27)
      g2 = ROOT.TF1('m2','gaus',27,34)
      g3 = ROOT.TF1('m3','gaus',32,40)
      total = ROOT.TF1('mstotal','gaus(0)+gaus(3)+gaus(6)+ \
                                  landau(9)+landau(12)+landau(15)',20,50)

      self.hist_container[name].Fit(g1,'R')
      self.hist_container[name].Fit(g2,'R+')
      self.hist_container[name].Fit(g3,'R+')

      par = [None]*18

      for i in range (0,3):
        par[i+0]  = g1.GetParameter(i)
        par[i+3]  = g2.GetParameter(i)
        par[i+6]  = g3.GetParameter(i)
        par[i+9]  = g4.GetParameter(i)
        par[i+12] = g5.GetParameter(i)
        par[i+15] = g6.GetParameter(i)

      print par

      self.hist_container[name].GetListOfFunctions().Clear()

      total.SetParameters(par[0], par[1], par[2],  \
                          par[3], par[4], par[5],  \
                          par[6], par[7], par[8],  \
                          par[9], par[10],par[11], \
                          par[12],par[13],par[14], \
                          par[15],par[16],par[17])
      fit=self.hist_container[name].Fit(total,'R+')

def print_data(data):
  pprint.pprint(data)

def Message(*args, **kwargs):
  string = ""
  for i in args:
    string = string + str(i)
  if "exc" in kwargs:
    if kwargs["exc"] == True or kwargs["exc"] == "testing":
      print string
  elif _config["verbose"] == True:
    print string
  log = open(_config["log_file"],"a")
  log.write(string)
  log.write("\n")
  log.close()
