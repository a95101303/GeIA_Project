#========Data-MC comparison for Electron=====#
import ROOT
import os
import math
from ROOT import TFile, TH1F, gDirectory, TCanvas, TPad, TProfile,TGraph, TGraphAsymmErrors,TGraphErrors
from ROOT import TH1D, TH1, TH1I
from ROOT import gStyle
from ROOT import gROOT
from ROOT import TStyle
from ROOT import TLegend
from ROOT import TMath
from ROOT import TPaveText
from ROOT import TLatex
from array import array

Hist_1=["min60m_ch0_hist","ped60m_ch0_hist","pedt60m_ch0_hist","diff_ped60m_ch0_pedt60m_ch0_hist"]
Hist_2=["min60m_ch0_ran_hist","ped60m_ch0_ran_hist","pedt60m_ch0_ran_hist","diff_ped60m_ch0_pedt60m_ch0_ran_hist"]
Hist_3=["min60m_ch0_af_hist","ped60m_ch0_af_hist","pedt60m_ch0_af_hist","diff_ped60m_ch0_pedt60m_ch0_af_hist"]

for i in range(0,4):
    f1= ROOT.TFile.Open("/Users/ms08962476/Report/AS_GeIA/Trailing_IoPAS/analysis/output_root/offset_cut_190627_2.root",'r')

    h1 = f1.Get(Hist_1[i])
    h2 = f1.Get(Hist_2[i])
    h3 = f1.Get(Hist_3[i])

    c = TCanvas("c1", "c1",0,0,500,500)
    gStyle.SetOptStat(0)
    leg = TLegend(0.4,0.1,0.8,0.2)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.SetBorderSize(0)
    leg.SetTextFont(22)
    leg.AddEntry(h1,"Random_trigger==0(Total)","p")
    leg.AddEntry(h2,"Random_trigger==1(OffSet)","p")


    if(i==0):
        h1.GetXaxis().SetRangeUser(-2000,33000)
        h1.GetYaxis().SetRangeUser(-300,-100)
        h1.SetXTitle("Max")
        h1.SetYTitle("Min")
        h1.SetZTitle("Arbitrary number")

    if(i==1):
        h1.GetXaxis().SetRangeUser(-2000,33000)
        #h1.GetYaxis().SetRangeUser(-221,-200)
        h1.GetYaxis().SetRangeUser(-210,-180)
        h1.SetXTitle("Max")
        h1.SetYTitle("Ped(First200bins)")
        h1.SetZTitle("Arbitrary number")
    if(i==2):
        h1.GetXaxis().SetRangeUser(-2000,33000)
        h1.GetYaxis().SetRangeUser(-210,-180)
        h1.SetXTitle("Max")
        h1.SetYTitle("Pedt(Last500bins)")
        h1.SetZTitle("Arbitrary number")
    if(i==3):
        h1.GetXaxis().SetRangeUser(-2000,33000)
        h1.GetYaxis().SetRangeUser(-20,15)
        h1.SetXTitle("Max")
        h1.SetYTitle("Ped-Pedt")
        h1.SetZTitle("Arbitrary number")

    h1.GetYaxis().SetLabelSize(0.03)
    h1.GetXaxis().SetTitleFont(22)
    h1.GetYaxis().SetTitleFont(22)
    h1.GetXaxis().SetLabelFont(22)
    h1.GetYaxis().SetLabelFont(22)

    h1.Draw("CLOZ")
    h2.Draw("CLOZsame")
    leg.Draw()

    c.Print(str(i)+"Try.pdf")







