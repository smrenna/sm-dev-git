// File: trimfile.cc
// Remove unwanted characters at the end of each line 
// and have exactly one empty line at the end of the file.
// Copyright (C) 2015 Torbjorn Sjostrand

// Stdlib header files.
#include <vector>
#include <cctype>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

// Used Stdlib elements.
using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

int main(int argc, char* argv[]) {

  // Check that correct number of command-line arguments
  if (argc != 2) {
    cout << " Wrong number of command-line arguments" << endl;
    return 1;
  }

  // Open input file.
  string fileName = argv[1];
  const char* fileCstring = fileName.c_str();
  ifstream is(fileCstring);  
  if (!is) {
    cout << " Input file " << fileName << " not found" << endl; 
    return 1;
  }

  // Read in input file.
  vector<std::string> lines;
  string line;
  while ( getline(is, line) ) lines.push_back(line);
  is.close();

  // Trim one line at a time.
  for (unsigned int i = 0; i < lines.size(); ++i) {
    line = lines[i];
    size_t posEnd = line.find_last_not_of(" \n\t\v\b\r\f\a");
    lines[i] = line.substr(0, posEnd + 1);
  }

  // Remove empty lines at the end.
  while (lines.back().length() == 0) lines.pop_back();

  // Open output file.
  ofstream os(fileCstring);  
  if (!os) {
    cout << " Output file " << fileName << " not found" << endl; 
    return 1;
  }

  // Write out all lines.
  for (unsigned int i = 0; i < lines.size(); ++i) os << lines[i] << endl;

  // Done.
  return 0;
}

