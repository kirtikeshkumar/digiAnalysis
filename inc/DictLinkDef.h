/*
** This file is required to create a root dictionary of the classes that have
** been defined.
** This allows for saving data in datatype of the classes
*/

#ifdef __ROOTCLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

// Add your classes here:
#pragma link C++ namespace digiAnalysis;

#pragma link C++ class digiAnalysis::WaveForm + ;
#pragma link C++ class digiAnalysis::singleHits + ;
#pragma link C++ class digiAnalysis::Pair + ;
#pragma link C++ class digiAnalysis::PSBar + ;
// add more as needed

#endif