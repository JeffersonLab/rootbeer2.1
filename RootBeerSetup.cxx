{
// Root macro to set up RootBeer, Bos, DST classes,
// and to set up for compilation of any functions using
// to sort Bos format data and rootDST format

//set include path for CINT to find inlcudes in macros
#pragma includepath ./ ./include


//set prompt
((TRint*)gROOT->GetApplication())->SetPrompt("rootbeer [%d] ");

// set the path for macros to be ./ and ./sample_code
gROOT->SetMacroPath(".:./sample_code:./include");

// load the mklib()  function which allows makes a shared library in
// $(TOP_DIR)/slib/Linux/libXX_C.so from <ANYDIR>/XX.C
gROOT->ProcessLine(".L MacroMaker.C");

// load any required libraries
//gSystem->Load("/usr/lib/libpthread.so");
gSystem->Load("$ROOTSYS/lib/libThread.so");
gSystem->Load("libRootBeer"); //Load the bos shared objects
gSystem->Load("libPhysics");
gSystem->Load("libeloss");
gSystem->Load("libExpTable");
//gSystem.Load("libPolHandler");
//~/slib/Linux64RHEL5/libDecayUtils.so 
gSystem->Load("~/packages/utilities/decayutils/depends/Linux64RHEL5/libDecayUtils.so");

//Set local includes in path

gSystem->SetIncludePath("-I./ -I./include");

}
