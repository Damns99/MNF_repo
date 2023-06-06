// File I/O for simple data in columns
#ifndef TEXT_IO_H
#define TEXT_IO_H

#include <fstream>
#include <string>
#include <cstdarg>
#include <sstream>
#include <vector>
#include <iomanip>

namespace textIo {
	
	namespace details {
		
		int getBetween(std::ifstream & infile, std::stringstream & sstr, char delim = '\t', char ignore = '#') {
			char ch;
			std::string temp = "";
			while (infile.get(ch)) {
				if (ch == ignore) std::getline(infile, temp);
				else if (ch == delim) return 1;
				else if (ch == '\n') return 1;
				else sstr << ch;
			}
			return 0;
		}
		
		template <typename T>
		int textInSlave(std::ifstream & infile, char delim, char ignore, unsigned int ind, T* dest);
		template <typename T, typename ... Args>
		int textInSlave(std::ifstream & infile, char delim, char ignore, unsigned int ind, T* dest, Args ... args);
		template <typename T>
		int textOutSlave(std::ofstream & outfile, char delim, unsigned int ind, T* orig);
		template <typename T, typename ... Args>
		int textOutSlave(std::ofstream & outfile, char delim, unsigned int ind, T* orig, Args ... args);
		
		int textInSlave(std::ifstream & infile, char delim, char ignore, unsigned int ind, std::vector<double>* dest);
		template <typename ... Args>
		int textInSlave(std::ifstream & infile, char delim, char ignore, unsigned int ind, std::vector<double>* dest, Args ... args);
		int textOutSlave(std::ofstream & outfile, char delim, unsigned int ind, std::vector<double>& orig);
		template <typename ... Args>
		int textOutSlave(std::ofstream & outfile, char delim, unsigned int ind, std::vector<double>& orig, Args ... args);
		int textOutSlave(std::ofstream & outfile, char delim, unsigned int ind, std::vector<std::vector<double>>& orig);
		template <typename ... Args>
		int textOutSlave(std::ofstream & outfile, char delim, unsigned int ind, std::vector<std::vector<double>>& orig, Args ... args);
		
		template <typename T>
		int textInSlave(std::ifstream & infile, char delim, char ignore, unsigned int ind, T* dest) {
			std::stringstream sstr;
			int endoffile = getBetween(infile, sstr, delim, ignore);
			sstr >> dest[ind];
			return endoffile;
		}
		template <typename T, typename ... Args>
		int textInSlave(std::ifstream & infile, char delim, char ignore, unsigned int ind, T* dest, Args ... args) {
			std::stringstream sstr;
			int endoffile = getBetween(infile, sstr, delim, ignore); 
			sstr >> dest[ind];
			if (endoffile == 0) return endoffile;
			else return textInSlave(infile, delim, ignore, ind, args ...);
		}
		int textInSlave(std::ifstream & infile, char delim, char ignore, unsigned int ind, std::vector<double>* dest) {
			std::stringstream sstr;
			int endoffile = getBetween(infile, sstr, delim, ignore);
			double tmp_;
			sstr >> tmp_;
			dest->push_back(tmp_);
			return endoffile;
		}
		template <typename ... Args>
		int textInSlave(std::ifstream & infile, char delim, char ignore, unsigned int ind, std::vector<double>* dest, Args ... args) {
			std::stringstream sstr;
			int endoffile = getBetween(infile, sstr, delim, ignore);
			double tmp_;
			sstr >> tmp_;
			dest->push_back(tmp_);
			if (endoffile == 0) return endoffile;
			else return textInSlave(infile, delim, ignore, ind, args ...);
		}
		
		template <typename T>
		int textOutSlave(std::ofstream & outfile, char delim, unsigned int ind, T* orig) {
			outfile << orig[ind] << '\n';
			return 0;
		}
		template <typename T, typename ... Args>
		int textOutSlave(std::ofstream & outfile, char delim, unsigned int ind, T* orig, Args ... args) {
			outfile << orig[ind] << ' ' << delim;
			return textOutSlave(outfile, delim, ind, args ...);
		}
		
		int textOutSlave(std::ofstream & outfile, char delim, unsigned int ind, std::vector<double>& orig) {
			outfile << std::setprecision(std::numeric_limits<double>::digits10);
			outfile << orig[ind] << '\n';
			return 0;
		}
		template <typename ... Args>
		int textOutSlave(std::ofstream & outfile, char delim, unsigned int ind, std::vector<double>& orig, Args ... args) {
			outfile << std::setprecision(std::numeric_limits<double>::digits10);
			outfile << orig[ind] << ' '  << delim;
			return textOutSlave(outfile, delim, ind, args ...);
		}
		int textOutSlave(std::ofstream & outfile, char delim, unsigned int ind, std::vector<std::vector<double>>& orig) {
			int i;
			for(i = 0; i < orig[0].size() - 1; i++) outfile << orig[ind][i] << ' '  << delim;
			outfile << orig[ind][i] << '\n';
			return 0;
		}
		template <typename ... Args>
		int textOutSlave(std::ofstream & outfile, char delim, unsigned int ind, std::vector<std::vector<double>>& orig, Args ... args) {
			for(auto& orig_i: orig[ind]) outfile << ' '  << orig_i;
			outfile << delim;
			return textOutSlave(outfile, delim, ind, args ...);
		}
	}
	
	template <typename ... Args>
	int textIn(const std::string filename, char delim = '\t', char ignore = '#', Args ... args) {
		std::ifstream infile;
		infile.open(filename, std::fstream::in);
		
		unsigned int ind = 0;
		while(details::textInSlave(infile, delim, ignore, ind, args ...) != 0) ind++;
		
		infile.close();
        return ind;
	}

	template <typename ... Args>
	void textOut(const std::string filename, char delim = '\t', char ignore = '#', const std::string header = "", unsigned int N = 0, bool append = false, Args ... args) {
		std::ofstream outfile;
		if (append) outfile.open(filename, std::fstream::out | std::fstream::app);
		else outfile.open(filename, std::fstream::out);
		
		outfile << ignore << header << '\n';
		
		for(unsigned int ind = 0; ind < N; ind++) details::textOutSlave(outfile, delim, ind, args ...);
		
		outfile.close();
	}
	
}

#endif