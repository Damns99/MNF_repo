// Command line argument parser
#ifndef CMDLINE_PARSER_H
#define CMDLINE_PARSER_H

#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <string.h>
#include <vector>
#include <variant>

#define CMDLINE_TYPES(_) _(char), _(short), _(int), _(long), _(unsigned), _(float), _(double), _(std::string)
#define CMDLINE_PARAMETER_TYPES(T) details::CmdlineParameter<T>
#define CMDLINE_POSSIBLE_TYPES CMDLINE_TYPES(CMDLINE_PARAMETER_TYPES)

#define HELP_RETURN 101

namespace cmdlineParser {
	
	namespace details {
		
		template <typename T>
		T parseCString(const char* cstring) {
			std::stringstream str;
			str << cstring;
			T res;
			str >> res;
			return res;
		}

		template <typename T>
		class CmdlineParameter{
			public:
				std::string name;
				T* valuePtr;
				std::string description;
				unsigned int nValues;
				
				CmdlineParameter(const std::string &new_name, T* new_pointer, const T &new_default, const std::string &new_description, unsigned int new_n_values) {
					name = new_name;
					valuePtr = new_pointer;
					nValues = new_n_values;
					for (unsigned int i = 0; i < new_n_values; i++) *(valuePtr + i) = new_default;
					description = new_description;
				}
		};
		
	}

	class CmdlineParser {
		private:
			std::vector<std::variant<CMDLINE_POSSIBLE_TYPES>> pos_;
			std::vector<std::variant<CMDLINE_POSSIBLE_TYPES>> opt_;
		
		public:
			template <typename T>
			void addPosParameter(const std::string &name, T* pointer, const T &default_value, const std::string &description) {
				details::CmdlineParameter<T> newparameter(name, pointer, default_value, description, 1);
				pos_.push_back(newparameter);
			}
			template <typename T>
			void addOptParameter(const std::string &name, T* pointer, const T &default_value, const std::string &description) {
				details::CmdlineParameter<T> newparameter(name, pointer, default_value, description, 1);
				opt_.push_back(newparameter);
			}
			template <typename T>
			void addPosArrayParameter(const std::string &name, T* pointer, const T &default_value, const std::string &description, unsigned int n_values) {
				details::CmdlineParameter<T> newparameter(name, pointer, default_value, description, n_values);
				pos_.push_back(newparameter);
			}
			template <typename T>
			void addOptArrayParameter(const std::string &name, T* pointer, const T &default_value, const std::string &description, unsigned int n_values) {
				details::CmdlineParameter<T> newparameter(name, pointer, default_value, description, n_values);
				opt_.push_back(newparameter);
			}
			
			void help(const char program_name[]) {
				std::cout << "Usage: " << program_name << " ";
				if (opt_.size() > 0) std::cout << "[options] ";
				for (auto& posi : pos_) std::visit([](auto& v){std::cout << v.name.c_str() << " ";}, posi);
				std::cout << "\nPositional:\n";
				for (auto& posi : pos_) std::visit([](auto& v){std::cout << std::left << std::setw(4) << "" << std::setw(16) << v.name.c_str() << std::setw(64) << v.description.c_str() << std::endl;}, posi);
				std::cout << "Options:\n";
				for (auto& opti : opt_) std::visit([](auto& v){std::cout << std::left << std::setw(4) << "" << std::setw(16) << v.name.c_str() << std::setw(64) << v.description.c_str() << std::endl;}, opti);
				std::cout << "\n";
			}
			
			void kickOff(const char program_name[]) {
				std::cout << "Starting " << program_name << " with arguments:" << std::endl;
				for (auto& posi : pos_) std::visit([](auto& v){
					std::cout << std::left << std::setw(4) << "" << std::setw(16) << v.name.c_str();
					std::cout << "[" << *(v.valuePtr);
					for (unsigned int k = 1; k < v.nValues; k++) std::cout << ", " << *(v.valuePtr + k);
					std::cout << "]" << std::endl;
				}, posi);
				for (auto& opti : opt_) std::visit([](auto& v){
					std::cout << std::left << std::setw(4) << "" << std::setw(16) << v.name.c_str();
					std::cout << "[" << *(v.valuePtr);
					for (unsigned int k = 1; k < v.nValues; k++) std::cout << ", " << *(v.valuePtr + k);
					std::cout << "]" << std::endl;
				}, opti);
				std::cout << "\n";
			}
			
			int parseAll(int argc, char* argv[]) {
				unsigned int pos_count = 0;
				for (int i = 1; i < argc; i++) {
					if (argv[i][0] == '-') {
						unsigned short found_ = false;
						for (auto& optj : opt_) std::visit([&](auto& v) mutable {
							if (v.name.compare(argv[i] + 1) == 0) {
								found_ = true;
								for (unsigned int k = 0; k < v.nValues && i + k + 1 < argc; k++) {
									std::stringstream str;
									str << argv[i + k + 1];
									str >> *(v.valuePtr + k);
								}
								i += v.nValues;
							}
						}, optj);
						if (!found_ && (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0)) {
							help(argv[0]);
							return HELP_RETURN;
						}
					}
					else {
						if (pos_count < pos_.size()) std::visit([&](auto& v) mutable {
							for (unsigned int k = 0; k < v.nValues && i + k < argc; k++) {
								std::stringstream str;
								str << argv[i + k];
								str >> *(v.valuePtr + k);
							}
							i += v.nValues - 1;
							pos_count++;
						}, pos_[pos_count]);
					}
				}
				return 0;
			}
	};
	
}

#endif
