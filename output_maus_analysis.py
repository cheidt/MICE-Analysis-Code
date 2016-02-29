import ROOT

class Output(object):
  def __init__(self, type):
    print "INITALIZING OUTPUTTER"
    print type
    self.hist_container = {}

  def Fill(self, name, title, *args, **kwargs):
    if "detector" in kwargs:
      title = title + " " + kwargs["detector"]
      name = name + "_" + kwargs["detector"]
    if "station" in kwargs:
      title = title + " " + str(kwargs["station"])
      name = name + "_" + str(kwargs["station"])
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
    out_root = ROOT.TFile("test_new_output.root",'RECREATE')
    for name in self.hist_container:
      self.hist_container[name].Write()
    out_root.Close()