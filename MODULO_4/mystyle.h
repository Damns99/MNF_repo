#ifndef MYSTYLE_H
#define MYSTYLE_H

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TROOT.h>

void myStyle() {
	TStyle* mystyle = new TStyle("mystyle", "mystyle");
	
	mystyle->SetLegendBorderSize(1);
	mystyle->SetLegendFillColor(kWhite);
	mystyle->SetFrameFillColor(kWhite);
	mystyle->SetCanvasColor(kWhite);
	mystyle->SetPadColor(kWhite);
	mystyle->SetFrameBorderMode(0);
	mystyle->SetCanvasBorderMode(0);
	mystyle->SetPadBorderMode(0);
	mystyle->SetPaperSize(20,26);
	
	mystyle->SetPadTopMargin(0.08);
	mystyle->SetPadRightMargin(0.05);
	mystyle->SetPadBottomMargin(0.1);
	mystyle->SetPadLeftMargin(0.1);
	
	mystyle->SetTextSize(0.03);
	mystyle->SetLabelSize(0.03, "xyz");
	mystyle->SetTitleSize(0.03, "xyz");
	mystyle->SetTitleFontSize(0.03);
	
	mystyle->SetTitleBorderSize(0);
	mystyle->SetTitleX(0.52);
	mystyle->SetTitleY(0.98);
	mystyle->SetTitleAlign(23);
	mystyle->SetTitleFillColor(kWhite);
	mystyle->SetTitleXOffset(1.4);
	mystyle->SetTitleYOffset(1.4);
	
	mystyle->SetMarkerStyle(8);
	mystyle->SetMarkerSize(0.4);
	mystyle->SetLineWidth(1);
	
	gROOT->SetStyle("mystyle");
	gROOT->ForceStyle();
}

#endif