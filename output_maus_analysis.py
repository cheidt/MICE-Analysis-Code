import ROOT

class Output:
  def Initialize (self):
    self.ROOT_Container = []

  def Fill(self, output):
    for key in output:
      print type(output[key])
      if type(output[key]) == "dict":
        pass