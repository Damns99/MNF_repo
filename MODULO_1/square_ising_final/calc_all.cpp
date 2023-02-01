#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <filesystem>
namespace fs = std::filesystem;

#include "bootstrap.h"
#include "text_io.h"
#include "cmdline_parser.h"

double mean(const std::vector<double>& vec) {
	double s = 0.;
	int length = vec.size();
	for(auto& ii: vec) s += ii / length;
	return s;
}

double derivative(const std::vector<double>& vec) {
	double s1 = 0., s2 = 0.;
	int length = vec.size();
	for(auto& ii: vec) {
		s1 += ii / length;
		s2 += (ii * ii) / length;
	}
	return s2 - s1 * s1;
}

double cumulBinder(const std::vector<double>& vec) {
	double s4 = 0., s2 = 0.;
	int length = vec.size();
	for(auto& ii: vec) {
		s4 += ii * ii * ii * ii / length;
		s2 += ii * ii / length;
	}
	return s4 / (s2 * s2);
}

double welfordVariance(const std::vector<double>& vec) {
	int count = 0;
	double mean = 0., M2 = 0.;
	for(auto& ii: vec) {
		count++;
		double delta = ii - mean;
		mean += delta / count;
		double delta2 = ii - mean;
		M2 += delta * delta2;
	}
	return M2 / count;
}

int main(int argc, char* argv[]) {
	long seed = -7745;
    std::string folder, outfilename, measfilename;
    std::vector<double> energy, magnetization, acceptance;
	int iscuda;
	
	cmdlineParser::CmdlineParser parser;
    parser.addPosParameter<std::string>("folder", &folder, "0220419144203", "[std::string] Working folder name.");
    parser.addPosParameter<std::string>("measfile", &measfilename, "metro_ising_meas.txt", "[std::string] Input file for the measures.");
    parser.addOptParameter<std::string>("outfile", &outfilename, "results.txt", "[std::string] Output file for the results.");
    parser.addOptParameter<int>("iscuda", &iscuda, 0, "[int] If != 0 searches in cudamesures/, if 0 in measures/.");
	if (parser.parseAll(argc, argv) == HELP_RETURN) return 0;
    parser.kickOff(argv[0]);
    
	if(iscuda != 0) fs::current_path(fs::current_path() / "cudameasures");
	else fs::current_path(fs::current_path() / "measures");
    fs::current_path(fs::current_path() / folder);
    int nmeas = textIo::textIn(measfilename, '\t', '#', &energy, &magnetization, &acceptance);
	for(auto& ii: magnetization) ii = abs(ii);
	
	std::string tmp;
	int lattice_length;
	std::ifstream tmpif;
	tmpif.open(measfilename, std::fstream::in);
	tmpif >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> lattice_length;
	tmpif.close();
    
    double energy_mean = mean(energy);
	double energy_error = bootstrap::bootstrapError(mean, energy, 256, 64, -std::time(NULL));
	double magnetization_mean = mean(magnetization);
	double magnetization_error = bootstrap::bootstrapError(mean, magnetization, 256, 64, -std::time(NULL));
    
    double energy_der = welfordVariance(energy) * lattice_length * lattice_length;
	double energy_der_error = bootstrap::bootstrapError(welfordVariance, energy, 256, 64, -std::time(NULL)) * lattice_length * lattice_length;
	double magnetization_der = welfordVariance(magnetization) * lattice_length * lattice_length;
	double magnetization_der_error = bootstrap::bootstrapError(welfordVariance, magnetization, 256, 64, -std::time(NULL)) * lattice_length * lattice_length;
	
	double energy_binder = cumulBinder(energy);
	double energy_binder_error = bootstrap::bootstrapError(cumulBinder, energy, 256, 64, -std::time(NULL));
	double magnetization_binder = cumulBinder(magnetization);
	double magnetization_binder_error = bootstrap::bootstrapError(cumulBinder, magnetization, 256, 64, -std::time(NULL));
	
	std::ofstream outfile;
	outfile.open(outfilename, std::fstream::out);
	outfile << "<E> = " << energy_mean << " +- " << energy_error << std::endl << "<|M|> = " << magnetization_mean << " +- " << magnetization_error << std::endl;
	outfile << "d<E>/dT = " << energy_der << " +- " << energy_der_error << std::endl << "d<|M|>/dh = " << magnetization_der << " +- " << magnetization_der_error << std::endl;
	outfile << "Binder(E) = " << energy_binder << " +- " << energy_binder_error << std::endl << "Binder(|M|) = " << magnetization_binder << " +- " << magnetization_binder_error << std::endl;
	outfile.close();
}