//functions to compile shared library or executables from a macro


void mklib(char *filename){
  mklib_with_flag(filename,0);
}
void mklib_deb(char *filename){
  mklib_with_flag(filename,1);
}
void mklib_with_flag(char *filename,int debug){
//turns <dir>/<file>.<suff> into $(TOPDIR)/slib/Linux/<file>_<suf>.so
//this makes a lib which can be loaded into root, and linked in the building of
//executables

  char *slibdir;
  char *slash;
  char *dot;
  int fnamelen;
  char suff[8];
  char filestub[100];
  char target[200];
  int result;
  char pwd[4]="./";
  char *libs;
    fprintf(stdout,"I am here");
  if(debug){
    char myMakeSharedLibString[500]="cd $BuildDir ; g++ -I$(TOP_DIR)/packages/utilities/decayutils -c $Opt -pipe -Wno-deprecated -Wno-unused -Wall -Woverloaded-virtual -fPIC -Iinclude -D_REENTRANT $IncludePath $SourceFiles ; g++ $Object
Files -shared -Wl,-soname,$LibName.so -g -o $SharedLib";
  }
  else{
    char myMakeSharedLibString[500]="cd $BuildDir ; g++ -I$(TOP_DIR)/packages/utilities/decayutils -c $Opt -pipe -Wno-deprecated -Wno-unused -Wall -Woverloaded-virtual -fPIC -Iinclude -D_REENTRANT $IncludePath $SourceFiles ; g++ $Object
Files -shared -Wl,-soname,$LibName.so -O -o $SharedLib";
  }

  char *makeSharedLibString;

  //get $ROOTBEER_SLIB, if not defined, make current dir
  if((slibdir = getenv("ROOTBEER_SLIB"))==NULL) slibdir=pwd;
  
  //find the last slash in the filename name or make start if no slash
  if((slash = strrchr(filename, '/'))!=NULL) slash++;
  else  slash = filename;

  //find the last dot in the filename and save _suffix, or bomb out
  if((dot = strrchr(filename, '.'))==NULL) return;
  sprintf(suff,"_%s",dot+1);
  
  fnamelen=(int)(dot-slash);                    //work out length from start to dot 
  sprintf(target,"%s/lib",slibdir);             //construct the target directory (just slibdir for now) 
  strncat(target,slash,fnamelen);               //append the filename
  //  strcat(target,suff);                      //append _suff
  slash = strrchr(target, '/');                 //get the new filename

  libs=gSystem->GetLibraries();
  if((strstr(libs,slash+1))!=NULL){             //see if loaded
    gSystem->Unload(slash+1);                   //unload if already loaded
  }
  makeSharedLibString=gSystem->GetMakeSharedLib();
  gSystem->SetMakeSharedLib(myMakeSharedLibString);

  //Now compile
  result= gSystem->CompileMacro(filename,"fk",target);
  // print some info for the user
  if(result==1){
    fprintf(stderr,"\nCONGRATULATIONS: Succussfully created the library as:\n%s.so\n\n",target);
    fprintf(stderr,"You can load it in BosSetup.cxx like this: gSystem.Load(\"lib%s\");\n",slash+4);
    fprintf(stderr,"And link an executable with it like this:  -L$(SODIR) -l%s\n\n",slash+4);    
  } 
  else{
    fprintf(stderr,"\nSORRY it didn't compile. Look carefully that the output above\n");
  }
  gSystem->SetMakeSharedLib(makeSharedLibString);
}

void mkexe(char *filename, char *extra){
  mkexe_with_flag(filename,extra,0);
}
void mkexe_deb(char *filename, char *extra){
  mkexe_with_flag(filename,extra,1);
}
void mkexe_with_flag(char *filename, char *extra, int debug){
//turns <dir>/<file>.<suff> into $(TOPDIR)/bin/Linux/<file>
//Makes an exe from a single macro file by running a system g++  command.
//eg. rootbeer>mkexe("sample_code/rbtest.C","-leloss")
//the "extra" arg can pass extra includes and libs
//-lRootReer and -lPhysics are already in.

  char *exedir;
  char *slash;
  char *dot;
  int fnamelen;
  char suff[8];
  char filestub[100];
  char target[200];
  int result;
  char pwd[4]="./";
  char *libs;

  char compileString[2500];

  //get $ROOTBEER_BIN, if not defined, make current dir
  if((exedir = getenv("ROOTBEER_BIN"))==NULL) exedir=pwd;
  
  //find the last slash in the filename name or make start if no slash
  if((slash = strrchr(filename, '/'))!=NULL) slash++;
  else  slash = filename;

  //find the last dot in the filename and save _suffix, or bomb out
  if((dot = strrchr(filename, '.'))==NULL) return;
  sprintf(suff,"_%s",dot+1);
  
  fnamelen=(int)(dot-slash);                    //work out length from start to dot 
  sprintf(target,"%s/",exedir);                 //construct the target directory (just exedir for now)
  strncat(target,slash,fnamelen);               //append the filename
  slash = strrchr(target, '/');                 //get the new filename
  
  if(debug){
  sprintf(compileString,"g++ -g -o %s -Wno-deprecated -Wall -Wno-unused -Woverloaded-virtual -fPIC %s -D_REENTRANT -DROOTEXE %s -lPhysics -lThread -lm -ldl -pthread -traditional -lRIO -lHist -lTree -L$ROOTBEER_SLIB -lRootBeer %s %s",tar
get, gSystem->GetIncludePath(),gSystem->GetLinkedLibs(),extra,filename);
  }
  else{
  sprintf(compileString,"g++ -O2 -o %s -Wno-deprecated -Wall -Wno-unused -Woverloaded-virtual -fPIC %s -D_REENTRANT -DROOTEXE %s -lPhysics -lThread -lm -ldl -pthread -traditional -lRIO -lHist -lTree -L$ROOTBEER_SLIB -lRootBeer %s %s",ta
rget, gSystem->GetIncludePath(),gSystem->GetLinkedLibs(),extra,filename);
  }

  fprintf(stdout,"%s\n",compileString);
  //Now compile
  result = system(compileString);
}

