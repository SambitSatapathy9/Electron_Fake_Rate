//DONT CHANGE THE FOLLOWING CODE AS THE PARAMETERS ARE SET
//CHECK

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooConstVar.h"
#include "RooBreitWigner.h"
#include "RooCrystalBall.h"
#include "TMath.h"
#include <iostream>
#include <TCanvas.h>

using namespace RooFit;

void ROOfitpTeg()
{
//Import input root file and read data
    TFile* infile = new TFile("testDEF.root", "READ");
    infile->cd();

//Define histogram and acces from the infile
    TH1F* histbwcb = (TH1F*)infile->Get("Z_ee200_220");
/*  TH1F* histbwcb = (TH1F*)infile->Get("Z_eg200_220");
    TH1F* histbwcb = (TH1F*)infile->Get("Z_ee220_270");
    TH1F* histbwcb = (TH1F*)infile->Get("Z_eg220_270");
    TH1F* histbwcb = (TH1F*)infile->Get("Z_ee270_320");
    TH1F* histbwcb = (TH1F*)infile->Get("Z_eg270_320");
    TH1F* histbwcb = (TH1F*)infile->Get("Z_ee320");
    TH1F* histbwcb = (TH1F*)infile->Get("Z_eg320");*/
    //Uncomment each line to access the input files

//Define observable "x"
    RooRealVar x("x","mass",60,120);

//Define RooDataHist - data from TH1F histogram is converted into a format suitable for statistical modeling using RooFit
    RooDataHist bwcbdata("bwcbdata","Data",x,Import(*histbwcb));

//Define nsig and nbkg so that the fit function and the data have comparable amount of events
    RooRealVar nsig("nsig","nsig",0,100000);
    RooRealVar nbkg("nbkg","nbkg",0,100000);

//BREIT WIGNER FUNCTION
    RooRealVar meanbw("meanbw","Mean of BW",91);
    RooRealVar widthbw("widthbw","Width of BW",2);

    RooBreitWigner breitwigner("breitwigner","Breit Wigner Function",x,meanbw,widthbw);

//CRYSTAL BALL FUNCTION
    RooRealVar meancb("meancb","Mean of CB",0.000075,-5,5);
    RooRealVar widthcb("widthcb","Width of CB",0.0001,-3,3);
    RooRealVar alpha("alpha","Alpha",0.0001,-10,10);
    RooRealVar n("n","Power law exponent",0.5,0,100);

    RooCrystalBall crystalball("crystalball","Crystal Ball Function",x,meancb,widthcb,alpha,n);

//Convolution of BW and CB Functions
    RooFFTConvPdf bwcb("bwcb","BWxCB",x,breitwigner,crystalball);

//CMS SHAPE (Sigmoid + Exponential)
    RooRealVar a("a","slope",0.0086);
    RooRealVar b("b","inflection point",0.0023);

    //RooGenericPdf arguments: name,title,formula,list of variables
    RooGenericPdf sigmoid("sigmoid","CMS SHAPE","((exp(b*x))/(1+exp((-a)*(x-91))))",RooArgSet(x,a,b));
//model final fit (BW*CB + CMS_SHAPE)
    RooAddPdf model("model","BWxCB + CMS_SHAPE",RooArgList(bwcb,sigmoid),RooArgList(nsig,nbkg));
    //It add two PDFs (bwcb,sigmoid) and (nsig,nbkg) so that the fit and data have events of comparable size

//Plots and Canvas
    gStyle->Print();
    gROOT->SetStyle("Plain");
    // gStyle->SetOptStat(1111);

    TCanvas* c1 = new TCanvas("c1","Convolution fit",1000,750);

    RooPlot* frame = x.frame();
//plot the data histogram to the frame
    bwcbdata.plotOn(frame,Name("bwcbdata"));

//fit model to data
    // model.fitTo(bwcbdata,Name("data"),Extended(),Range(60,120));
    model.fitTo(bwcbdata,Name("data"),Range(60,120));

//plot the model and its components to frame
    model.plotOn(frame,Name("model"),LineColor(kRed));
    model.plotOn(frame,Name("sigmoid"),Components(sigmoid),LineStyle(kDashed));
    model.plotOn(frame,Name("bwcb"),Components(bwcb),LineColor(kGreen));

//set frame title
    frame->SetTitle("Z(ee)InvariantMass200_220");
 /* frame->SetTitle("Z(eg)InvariantMass200_220");
    frame->SetTitle("Z(ee)InvariantMass220_270");
    frame->SetTitle("Z(eg)InvariantMass220_270");
    frame->SetTitle("Z(ee)InvariantMass270_320");
    frame->SetTitle("Z(eg)InvariantMass270_320");
    frame->SetTitle("Z(ee)InvariantMass320");
    frame->SetTitle("Z(eg)InvariantMass320");*/
    //Uncomment each line to set the title of the plots
    frame->GetXaxis()->SetTitle("Mass in GeV");
    frame->GetYaxis()->SetTitle("Number of events");
    frame->Draw();

//Adding Legends to the plots
    TLegend* leg = new TLegend(0.60,0.45,0.88,0.87);
    leg->SetFillColor(kWhite);
    leg->SetLineColor(kWhite);
    leg->AddEntry("model","BWxCB + CMS_SHAPE");
    leg->AddEntry("sigmoid","Cms Shape(Sigmoid+Exp)");
    leg->AddEntry("bwcb","BWxCB");


//////////////************************************/////////////////////
    // After the model is plotted on the frame:

    // Get the parameter values
    double alphaVal = alpha.getVal();
    double meanbwVal = meanbw.getVal();
    double widthbwVal = widthbw.getVal();
    double meancbVal = meancb.getVal();
    double widthcbVal = widthcb.getVal();
    double nsigVal = nsig.getVal();
    double nbkgVal = nbkg.getVal();
    double aVal = a.getVal();
    double bVal = b.getVal();
    double nVal = n.getVal();
    // Convert the parameter values to strings
    std::string alphaStr = Form("#alpha_{CB} = %.3f", alphaVal);
    std::string meancbStr = Form("meanCB = %.3f", meancbVal);
    std::string widthcbStr = Form("widthCB = %.3f", widthcbVal);
    std::string nStr = Form("n_{CB} = %0.3f",nVal);

    std::string meanbwStr = Form("meanBW = %.3f", meanbwVal);
    std::string widthbwStr = Form("widthBW = %.3f", widthbwVal);

    std::string slopeStr = Form("slope = %.5f",aVal);
    std::string inflectionStr = Form("inflection = %.5f",bVal);

    std::string nsigStr = Form("nsig = %.0f",nsigVal);
    std::string nbkgStr = Form("nbkg = %.0f", nbkgVal);

    // Add the parameter values to the legend
    // leg->AddEntry((TObject*)0, meanbwStr.c_str(), "");
    // leg->AddEntry((TObject*)0, widthbwStr.c_str(), "");

    leg->AddEntry((TObject*)0, alphaStr.c_str(), "");
    leg->AddEntry((TObject*)0, meancbStr.c_str(), "");
    leg->AddEntry((TObject*)0, widthcbStr.c_str(), "");
    leg->AddEntry((TObject*)0, nStr.c_str(),"");

    leg->AddEntry((TObject*)0,slopeStr.c_str(),"");
    leg->AddEntry((TObject*)0,inflectionStr.c_str(),"");

    leg->AddEntry((TObject*)0, nsigStr.c_str(),"");
    leg->AddEntry((TObject*)0, nbkgStr.c_str(), "");

//////////////************************************/////////////////////

    leg->Draw();

    c1->Draw();
    c1->Update();
    c1->SaveAs("Zee200_220.png");
/*  c1->SaveAs("Zeg200_220.png");
    c1->SaveAs("Zee220_270.png");
    c1->SaveAs("Zeg220_270.png");
    c1->SaveAs("Zee270_320.png");
    c1->SaveAs("Zeg270_320.png");
    c1->SaveAs("Zee320.png");
    c1->SaveAs("Zeg320.png");*/
    //uncomment each line above to save the plots
}
