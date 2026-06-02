// A simple event counter from the Condor report files
// This allows one to quickly estimate the number of events
// sent by the STAR scheduler to be processed, as can be 
// compared against the number of events analyzed by code/
// known to be contained in a data sample. 
//
// To compare against the STAR database/records, visit the run log:
//
//  (All DAQ events) 
//	https://online.star.bnl.gov/RunLog/
//  (production release events)
//      https://www.star.bnl.gov/public/comp/prod/ProdSumPicoDst.html
// Written by Ryan Hamilton on 4/20/2026

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

int count_events(const char* job_identifier) {
  bool _debug = false;
  
  char filename[150];
  snprintf(filename, 150, "report/sched%s.report", job_identifier);
  
  int count_table_delims = 0;
  std::string table_delim = "+----------------";

  std::ifstream infile(filename);
  std::string line;
  int total_events, total_files, total_processes = 0;
  while (getline(infile, line)) {
    if (table_delim.compare(line.substr(0, table_delim.size())) == 0) {
      ++count_table_delims;
      if (count_table_delims <= 3) continue;
      else break;
    } else if (count_table_delims < 3) continue;
    
    std::stringstream linestream(line);
    std::string entry;
    for (int i = 0; i < 3; ++i) {
      std::getline(linestream, entry, '|');
      if (_debug) std::cout << "Skipping string [" << entry << "]." << std::endl;
    }
    
    // Add number of files
    std::getline(linestream, entry, '|');
    if (_debug) std::cout << "attempting to convert [" << entry << "]." << std::endl;
    total_files += std::stoi(entry);
    
    // Add number of events
    std::getline(linestream, entry, '|');
    if (_debug) std::cout << "attempting to convert [" << entry << "]." << std::endl;
    total_events += std::stoi(entry);
    
    ++total_processes;
  }

  std::cout << "Processed report for job ID :: " << job_identifier << std::endl;
  std::cout << "   -> \033[31m" << total_processes << "\033[39m processes." << std::endl;
  std::cout << "   -> \033[31m" << total_files << "\033[39m files." << std::endl;
  std::cout << "   -> \033[31m" << total_events << "\033[39m events." << std::endl;
  return total_events;
}

int main(int argc, char *argv[] ) {
  if (argc <= 1) {
    std::cerr << "Insufficient arguments!" << std::endl;
    std::cout << "arg1 :: job identifier to use for analyzing the report file." << std::endl;
    return 1;
  }
  
  count_events(argv[1]);

  return 0;
}
