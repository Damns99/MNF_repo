#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include <filesystem>
namespace fs = std::filesystem;

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLegend.h>

#include "mystyle.h"

#include "cmdline_parser.h"
#include "integrate_LIF.h"
#include "text_io.h"

int main(int argc, char* argv[]) {
	myStyle();
	
	std::string filenamelist[5];
	double miny = 0.86, maxy = 1.02;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosArrayParameter<std::string>("filenamelist", filenamelist, "", "output file name list with extensions [string[5]]", 5);
	
	if(parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
	
	TCanvas* canvas = new TCanvas("canvas0", "Canvas0", 600, 400);
	
	TMultiGraph* multigraph = new TMultiGraph();
	std::string title = "";
	
	fs::current_path(fs::current_path() / "measures");
	
	for(int ff = 0; ff < 5; ff++) {
		std::string filename = filenamelist[ff];
		if(filename == "") break;
		std::ifstream infile;
		infile.open(filename, std::fstream::in);
		
		std::string tmpline = "";
		char ch;
		infile.get(ch); std::getline(infile, title);
		std::getline(infile, tmpline);
		
		while(infile.peek() != EOF) {
			
			std::vector<double> x, y;
			
			for(int ii = 0; ii < 1e3 && infile.peek() != EOF; ii++) {
				double xx, yy;
				infile >> xx >> yy;
				x.push_back(xx); y.push_back(yy);
			}
			
			TGraph* graph = new TGraph(x.size(), x.data(), y.data());
			graph->SetMarkerColor(ff+1);
			graph->SetLineColor(ff+1);
			multigraph->Add((TGraph*) graph->Clone(), "");
		}
		
		infile.close();
	}
	
	canvas->cd();
	canvas->SetGrid();
	multigraph->SetTitle(";t [Cm/g];V [Vrest]");
	multigraph->SetMaximum(maxy);
	multigraph->SetMinimum(miny);
	multigraph->Draw("APL");
	canvas->SaveAs("test_integration_comparison.pdf");
	
}