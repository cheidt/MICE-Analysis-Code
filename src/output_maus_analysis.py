import ROOT
import os
from config_maus_analysis import out_config as _config

class Output(object):
  def __init__(self, type):
    self.name = type
    self.hist_container = {}

  def Fill(self, name, title, *args, **kwargs):
    if "detector" in kwargs:
      title = title + " " + kwargs["detector"]
      name  = name + "_" + kwargs["detector"]
    if "station" in kwargs:
      title = title + " " + str(kwargs["station"])
      name  = name + "_" + str(kwargs["station"])
    if "iteration" in kwargs:
      title = title + " itr " + str(kwargs["iteration"])
      name  = name + "_" + str(kwargs["iteration"])
    if name not in self.hist_container:
      self.Initialize(name, title, args)
    else:
      if len(args) == 8:
        self.hist_container[name].Fill(args[0], args[1])
      else:
        self.hist_container[name].Fill(args[0])

  def Initialize(self, name, title, args):
    if len(args) == 8:
      self.hist_container[name] = ROOT.TH2D(name, title, args[2], args[3], \
                                                         args[4], args[5], \
                                                         args[6], args[7])
    else:
      self.hist_container[name] = ROOT.TH1D(name, title,args[1], args[2], \
                                                        args[3])

  def Write(self):
    file_name = _config["output_dir"] + self.name + _config["output_file"]
    print "WRITTING ANALYSIS FILE: ",file_name
    out_root = ROOT.TFile(file_name,'RECREATE')
    for name in self.hist_container:
      self.hist_container[name].Write()
    out_root.Close()

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
  