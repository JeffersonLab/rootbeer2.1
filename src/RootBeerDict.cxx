// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME srcdIRootBeerDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "include/TRootBeer.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_TRootBeer(void *p = nullptr);
   static void *newArray_TRootBeer(Long_t size, void *p);
   static void delete_TRootBeer(void *p);
   static void deleteArray_TRootBeer(void *p);
   static void destruct_TRootBeer(void *p);
   static void streamer_TRootBeer(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRootBeer*)
   {
      ::TRootBeer *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRootBeer >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TRootBeer", ::TRootBeer::Class_Version(), "include/TRootBeer.h", 24,
                  typeid(::TRootBeer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TRootBeer::Dictionary, isa_proxy, 16,
                  sizeof(::TRootBeer) );
      instance.SetNew(&new_TRootBeer);
      instance.SetNewArray(&newArray_TRootBeer);
      instance.SetDelete(&delete_TRootBeer);
      instance.SetDeleteArray(&deleteArray_TRootBeer);
      instance.SetDestructor(&destruct_TRootBeer);
      instance.SetStreamerFunc(&streamer_TRootBeer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRootBeer*)
   {
      return GenerateInitInstanceLocal(static_cast<::TRootBeer*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TRootBeer*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TRootBeer::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TRootBeer::Class_Name()
{
   return "TRootBeer";
}

//______________________________________________________________________________
const char *TRootBeer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootBeer*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TRootBeer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRootBeer*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TRootBeer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootBeer*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TRootBeer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRootBeer*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TRootBeer::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRootBeer.

   TObject::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRootBeer(void *p) {
      return  p ? new(p) ::TRootBeer : new ::TRootBeer;
   }
   static void *newArray_TRootBeer(Long_t nElements, void *p) {
      return p ? new(p) ::TRootBeer[nElements] : new ::TRootBeer[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRootBeer(void *p) {
      delete (static_cast<::TRootBeer*>(p));
   }
   static void deleteArray_TRootBeer(void *p) {
      delete [] (static_cast<::TRootBeer*>(p));
   }
   static void destruct_TRootBeer(void *p) {
      typedef ::TRootBeer current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TRootBeer(TBuffer &buf, void *obj) {
      ((::TRootBeer*)obj)->::TRootBeer::Streamer(buf);
   }
} // end of namespace ROOT for class ::TRootBeer

namespace {
  void TriggerDictionaryInitialization_RootBeerDict_Impl() {
    static const char* headers[] = {
"include/TRootBeer.h",
nullptr
    };
    static const char* includePaths[] = {
"/cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/almalinux9-gcc11/local/root/6.30.04/include/",
"/u/home/nicholas/HyperonDeut/rootbeer2.1/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "RootBeerDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$include/TRootBeer.h")))  TRootBeer;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "RootBeerDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "include/TRootBeer.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"TRootBeer", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("RootBeerDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_RootBeerDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_RootBeerDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_RootBeerDict() {
  TriggerDictionaryInitialization_RootBeerDict_Impl();
}
