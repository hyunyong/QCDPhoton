from ROOT import *
gROOT.SetBatch(1)
gStyle.SetOptStat(0)

tf = TFile("QCDPhoton.root")
tr = tf.Get("QCDPhoton")

outFile = open("photonM.txt", "w")
for e in tr:
  print "Event number: ", e.evn
  print "  Number of jet: ", e.nJet
  for i, pJet in enumerate(e.photon):
     print "    Number of photon at jet %d: "%i, len(pJet)
     for j, p in enumerate(pJet):
       print "      %d photon 4 vec.: "%j, p.Px(), p.Py(), p.Pz(), p.Pt()
       outFile.write("%f %f %f %f\n"%(p.Px(), p.Py(), p.Pz(), p.E()))
