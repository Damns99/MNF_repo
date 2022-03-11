// File I/O for simple data in columns
#ifndef TEXT_IO_H
#define TEXT_IO_H

#include <fstream>
#include <string>
#include <cstdarg>
#include <sstream>

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
		
		template <typename T>
		int textOutSlave(std::ofstream & outfile, char delim, unsigned int ind, T* orig) {
			outfile << orig[ind] << '\n';
			return 0;
		}
		template <typename T, typename ... Args>
		int textOutSlave(std::ofstream & outfile, char delim, unsigned int ind, T* orig, Args ... args) {
			outfile << orig[ind] << delim;
			return textOutSlave(outfile, delim, ind, args ...);
		}
		
	}
	
	template <typename ... Args>
	void textIn(const std::string filename, char delim = '\t', char ignore = '#', Args ... args) {
		std::ifstream infile;
		infile.open(filename, std::fstream::in);
		
		unsigned int ind = 0;
		while(details::textInSlave(infile, delim, ignore, ind, args ...) != 0) ind++;
		
		infile.close();
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