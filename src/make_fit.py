#!/usr/bin/env python
import ROOT

class Fit():
  print "Loading previously used ROOT file"
  fit_file = ROOT.TFile("~/Dropbox/Work/Scripts/output/Mass_output_140_Diff0_lattice1_5_LiH.root", "READ")
  # plot = fit_file.Get("Mass_upstream")
  plot = fit_file.Get("Mass_downstream")

  fit_file.Draw()

  # raw_input("Press Enter to Exit")

  g1 = ROOT.TF1('m1', 'gaus', 80, 130)
  g2 = ROOT.TF1('m2', 'landau', 90, 130)
  g3 = ROOT.TF1('m3', 'gaus', 120, 180)
  # g4 = ROOT.TF1('m4', 'landau', 130, 160)
  # g5 = ROOT.TF1('m2', 'landau', 27, 34)
  # g6 = ROOT.TF1('m3', 'landau', 32, 40)
  # total = ROOT.TF1('mstotal','gaus(0)+landau(3)',85,130)
  total = ROOT.TF1('mstotal','gaus(0)+landau(3)+gaus(6)',80,180)

  # g3 = ROOT.TF1('m3', 'gaus', 32, 40)
  # total = ROOT.TF1('mstotal', 'gaus(0)+gaus(3)+gaus(6)+ \
  #                              landau(9)+landau(12)+landau(15)', 20, 50)

  plot.Fit(g1, 'R')
  plot.Fit(g2, 'R+')
  plot.Fit(g3, 'R+')
  # plot.Fit(g4, 'R+')

  par = [None] * 9

  for i in range(0, 3):
    par[i + 0] = g1.GetParameter(i)
    par[i + 3] = g2.GetParameter(i)
    par[i + 6] = g3.GetParameter(i)
    # par[i + 9] = g4.GetParameter(i)


  print par

  plot.GetListOfFunctions().Clear()

  total.SetParameters(par[0], par[1], par[2], \
                      par[3], par[4], par[5], \
                      par[6], par[7], par[8])

  # total.SetParameters(par[0], par[1], par[2], \
  #                     par[3], par[4], par[5])
  fit = plot.Fit(total, 'R+')

  tot_par = [None] * 9
  for i in range(len(tot_par)):
    tot_par[i] = total.GetParameter(i)

  print ROOT.Math.landau_pdf(tot_par[3], tot_par[4], tot_par[5])

  raw_input("Press Enter to Exit")

if __name__ == "__main__":
  Fit()