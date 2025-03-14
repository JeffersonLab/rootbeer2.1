//
// File generated by /apps/root/5.34.36/root//bin/rootcint at Mon Nov 25 15:22:46 2019

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME srcdIBosDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "BosDict.h"

#include "TClass.h"
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

// Direct notice to TROOT of the dictionary's loading.
namespace {
   static struct DictInit {
      DictInit() {
         ROOT::RegisterModule();
      }
   } __TheDictionaryInitializer;
}

// START OF SHADOWS

namespace ROOTShadow {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOTShadow
// END OF SHADOWS

namespace ROOTDict {
   void TBos_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void delete_TBos(void *p);
   static void deleteArray_TBos(void *p);
   static void destruct_TBos(void *p);
   static void streamer_TBos(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static ROOT::TGenericClassInfo *GenerateInitInstanceLocal(const ::TBos*)
   {
      ::TBos *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TBos >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TBos", ::TBos::Class_Version(), "./include/TBos.h", 16,
                  typeid(::TBos), ::ROOT::DefineBehavior(ptr, ptr),
                  &::TBos::Dictionary, isa_proxy, 0,
                  sizeof(::TBos) );
      instance.SetDelete(&delete_TBos);
      instance.SetDeleteArray(&deleteArray_TBos);
      instance.SetDestructor(&destruct_TBos);
      instance.SetStreamerFunc(&streamer_TBos);
      return &instance;
   }
   ROOT::TGenericClassInfo *GenerateInitInstance(const ::TBos*)
   {
      return GenerateInitInstanceLocal((::TBos*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TBos*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOTDict

//______________________________________________________________________________
atomic_TClass_ptr TBos::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TBos::Class_Name()
{
   return "TBos";
}

//______________________________________________________________________________
const char *TBos::ImplFileName()
{
   return ::ROOTDict::GenerateInitInstanceLocal((const ::TBos*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TBos::ImplFileLine()
{
   return ::ROOTDict::GenerateInitInstanceLocal((const ::TBos*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void TBos::Dictionary()
{
   fgIsA = ::ROOTDict::GenerateInitInstanceLocal((const ::TBos*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *TBos::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gCINTMutex); if(!fgIsA) {fgIsA = ::ROOTDict::GenerateInitInstanceLocal((const ::TBos*)0x0)->GetClass();} }
   return fgIsA;
}

//______________________________________________________________________________
void TBos::Streamer(TBuffer &R__b)
{
   // Stream an object of class TBos.

   TRootBeer::Streamer(R__b);
}

//______________________________________________________________________________
void TBos::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class TBos.
      TClass *R__cl = ::TBos::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fDataThread", &fDataThread);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fScannerThread", &fScannerThread);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fERing", &fERing);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fBankList[220][5]", fBankList);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fBankType[220]", fBankType);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fDataFile[80]", fDataFile);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fNanalysed", &fNanalysed);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fSwapFlag", &fSwapFlag);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fRingEvent", &fRingEvent);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fPrevRingEvent", &fPrevRingEvent);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fServedBankNo", &fServedBankNo);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fBufferStatus[2]", fBufferStatus);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fBufferNread[2]", fBufferNread);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fRBuff[2]", &fRBuff);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fKillScanner", &fKillScanner);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fKillDataServer", &fKillDataServer);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fIsScanner", &fIsScanner);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fIsDataServer", &fIsDataServer);
      TRootBeer::ShowMembers(R__insp);
}

namespace ROOTDict {
   // Wrapper around operator delete
   static void delete_TBos(void *p) {
      delete ((::TBos*)p);
   }
   static void deleteArray_TBos(void *p) {
      delete [] ((::TBos*)p);
   }
   static void destruct_TBos(void *p) {
      typedef ::TBos current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TBos(TBuffer &buf, void *obj) {
      ((::TBos*)obj)->::TBos::Streamer(buf);
   }
} // end of namespace ROOTDict for class ::TBos

/********************************************************
* src/BosDict.cxx
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableBosDict();

extern "C" void G__set_cpp_environmentBosDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("include/TBos.h");
  G__cpp_reset_tagtableBosDict();
}
#include <new>
extern "C" int G__cpp_dllrevBosDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* TBos */
static int G__BosDict_255_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TBos* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 3
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new TBos(
(const char*) G__int(libp->para[0]), (void*) G__int(libp->para[1])
, (Int_t) G__int(libp->para[2]));
   } else {
     p = new((void*) gvp) TBos(
(const char*) G__int(libp->para[0]), (void*) G__int(libp->para[1])
, (Int_t) G__int(libp->para[2]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__BosDictLN_TBos));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TBos::WordSwap((char*) G__int(libp->para[0]), (int) G__int(libp->para[1]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TBos::DataSwap((char*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 89, (long) TBos::DataServer((void*) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 89, (long) TBos::BankScanner((void*) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TBos::printRHead((RecHead_t*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TBos::printRSHead((recSegHead_t*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TBos::printDataSHead((dataSegHead_t*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) TBos::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TBos::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) TBos::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TBos::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TBos*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TBos::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TBos::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TBos::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__BosDict_255_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TBos::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__BosDict_255_0_23(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   TBos* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new TBos(*(TBos*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__BosDictLN_TBos));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef TBos G__TTBos;
static int G__BosDict_255_0_24(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (TBos*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((TBos*) (soff+(sizeof(TBos)*i)))->~G__TTBos();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (TBos*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((TBos*) (soff))->~G__TTBos();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__BosDict_255_0_25(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TBos* dest = (TBos*) G__getstructoffset();
   *dest = *(TBos*) libp->para[0].ref;
   const TBos& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* TBos */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncBosDict {
 public:
  G__Sizep2memfuncBosDict(): p(&G__Sizep2memfuncBosDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncBosDict::*p)();
};

size_t G__get_sizep2memfuncBosDict()
{
  G__Sizep2memfuncBosDict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceBosDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__BosDictLN_TBos))) {
     TBos *G__Lderived;
     G__Lderived=(TBos*)0x1000;
     {
       TRootBeer *G__Lpbase=(TRootBeer*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__BosDictLN_TBos),G__get_linked_tagnum(&G__BosDictLN_TRootBeer),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__BosDictLN_TBos),G__get_linked_tagnum(&G__BosDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableBosDict() {

   /* Setting up typedef entry */
   G__search_typename2("Int_t",105,-1,0,-1);
   G__setnewtype(-1,"Signed integer 4 bytes (int)",0);
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__BosDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BosDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BosDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BosDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BosDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__BosDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BosDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BosDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BosDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BosDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<UInt_t>",117,G__get_linked_tagnum(&G__BosDictLN_vectorlEunsignedsPintcOallocatorlEunsignedsPintgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TString>",117,G__get_linked_tagnum(&G__BosDictLN_vectorlETStringcOallocatorlETStringgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__BosDictLN_reverse_iteratorlEvectorlETStringcOallocatorlETStringgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BosDictLN_vectorlETStringcOallocatorlETStringgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__BosDictLN_reverse_iteratorlEvectorlETStringcOallocatorlETStringgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__BosDictLN_vectorlETStringcOallocatorlETStringgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* TBos */
static void G__setup_memvarTBos(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__BosDictLN_TBos));
   { TBos *p; p=(TBos*)0x1000; if (p) { }
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__BosDictLN_TThread),-1,-1,4,"fDataThread=",0,"thread pointer");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__BosDictLN_TThread),-1,-1,4,"fScannerThread=",0,"thread pointer");
   G__memvar_setup((void*)((long)(&p->fERing)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__BosDictLN_eventInf_t),-1,-1,1,"fERing=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->fBankList)-(long)(p)),99,0,0,-1,G__defined_typename("Char_t"),-1,1,"fBankList[220][5]=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->fBankType)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"fBankType[220]=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->fDataFile)-(long)(p)),99,0,0,-1,G__defined_typename("Char_t"),-1,1,"fDataFile[80]=",0,"hold the name of the current data file ");
   G__memvar_setup((void*)((long)(&p->fNanalysed)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"fNanalysed=",0,"count no of analysed events");
   G__memvar_setup((void*)((long)(&p->fSwapFlag)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"fSwapFlag=",0,"flag to indicate byte swapping");
   G__memvar_setup((void*)((long)(&p->fRingEvent)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"fRingEvent=",0,"current ring event index");
   G__memvar_setup((void*)((long)(&p->fPrevRingEvent)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"fPrevRingEvent=",0,"previous ring event index");
   G__memvar_setup((void*)((long)(&p->fServedBankNo)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"fServedBankNo=",0,"count the no of banks to be served");
   G__memvar_setup((void*)((long)(&p->fBufferStatus)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"fBufferStatus[2]=",0,"is the buffer READ, SCANNED, or ANALYSED ");
   G__memvar_setup((void*)((long)(&p->fBufferNread)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"fBufferNread[2]=",0,"no of bytes read into buffer ");
   G__memvar_setup((void*)((long)(&p->fRBuff)-(long)(p)),67,0,0,-1,-1,-1,1,"fRBuff[2]=",0,"Two buffers to take records read");
   G__memvar_setup((void*)((long)(&p->fKillScanner)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"fKillScanner=",0,"Flag to end the theads");
   G__memvar_setup((void*)((long)(&p->fKillDataServer)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"fKillDataServer=",0,"Flag to end the theads");
   G__memvar_setup((void*)((long)(&p->fIsScanner)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"fIsScanner=",0,"Flag to say if thread running");
   G__memvar_setup((void*)((long)(&p->fIsDataServer)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"fIsDataServer=",0,"ditto");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__BosDictLN_TClass),G__defined_typename("atomic_TClass_ptr"),-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarBosDict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncTBos(void) {
   /* TBos */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__BosDictLN_TBos));
   G__memfunc_setup("TBos",376,G__BosDict_255_0_1, 105, G__get_linked_tagnum(&G__BosDictLN_TBos), -1, 0, 3, 1, 1, 0, 
"C - - 10 - - Y - - 0 - - "
"i - 'Int_t' 0 - -", "class constructor", (void*) NULL, 0);
   G__memfunc_setup("GetEvent",802,(G__InterfaceMethod) NULL,105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "main function called for each event", (void*) NULL, 1);
   G__memfunc_setup("WordSwap",823,G__BosDict_255_0_3, 121, -1, -1, 0, 2, 3, 1, 0, 
"C - - 0 - - i - - 0 - -", "if written on big(little)endian and read on little(big)...", (void*) G__func2void( (void (*)(char*, int))(&TBos::WordSwap) ), 0);
   G__memfunc_setup("DataSwap",789,G__BosDict_255_0_4, 121, -1, -1, 0, 1, 3, 1, 0, "C - - 0 - -", "...endian, some byte swapping needs to be done", (void*) G__func2void( (void (*)(char*))(&TBos::DataSwap) ), 0);
   G__memfunc_setup("StartServer",1157,(G__InterfaceMethod) NULL,105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "initialise some things after constructor called", (void*) NULL, 1);
   G__memfunc_setup("DataServer",1009,G__BosDict_255_0_6, 89, -1, -1, 0, 1, 3, 1, 0, "Y - - 0 - arg", "some debugging print functions for bos files", (void*) G__func2void( (void* (*)(void*))(&TBos::DataServer) ), 0);
   G__memfunc_setup("BankScanner",1094,G__BosDict_255_0_7, 89, -1, -1, 0, 1, 3, 1, 0, "Y - - 0 - arg", (char*)NULL, (void*) G__func2void( (void* (*)(void*))(&TBos::BankScanner) ), 0);
   G__memfunc_setup("printRHead",1009,G__BosDict_255_0_8, 121, -1, -1, 0, 1, 3, 1, 0, "U 'RecHead_t' - 0 - s", (char*)NULL, (void*) G__func2void( (void (*)(RecHead_t*))(&TBos::printRHead) ), 0);
   G__memfunc_setup("printRSHead",1092,G__BosDict_255_0_9, 121, -1, -1, 0, 1, 3, 1, 0, "U 'recSegHead_t' - 0 - s", (char*)NULL, (void*) G__func2void( (void (*)(recSegHead_t*))(&TBos::printRSHead) ), 0);
   G__memfunc_setup("printDataSHead",1388,G__BosDict_255_0_10, 121, -1, -1, 0, 1, 3, 1, 0, "U 'dataSegHead_t' - 0 - s", (char*)NULL, (void*) G__func2void( (void (*)(dataSegHead_t*))(&TBos::printDataSHead) ), 0);
   G__memfunc_setup("Class",502,G__BosDict_255_0_11, 85, G__get_linked_tagnum(&G__BosDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&TBos::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__BosDict_255_0_12, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TBos::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__BosDict_255_0_13, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&TBos::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__BosDict_255_0_14, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&TBos::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__BosDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__BosDict_255_0_18, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__BosDict_255_0_19, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TBos::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__BosDict_255_0_20, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TBos::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__BosDict_255_0_21, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TBos::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__BosDict_255_0_22, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TBos::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("TBos", 376, G__BosDict_255_0_23, (int) ('i'), G__get_linked_tagnum(&G__BosDictLN_TBos), -1, 0, 1, 1, 1, 0, "u 'TBos' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~TBos", 502, G__BosDict_255_0_24, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__BosDict_255_0_25, (int) ('u'), G__get_linked_tagnum(&G__BosDictLN_TBos), -1, 1, 1, 1, 1, 0, "u 'TBos' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncBosDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {
}

static void G__cpp_setup_global2() {
}

static void G__cpp_setup_global3() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalBosDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
  G__cpp_setup_global2();
  G__cpp_setup_global3();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {
}

static void G__cpp_setup_func13() {
}

static void G__cpp_setup_func14() {
}

static void G__cpp_setup_func15() {
}

static void G__cpp_setup_func16() {
}

static void G__cpp_setup_func17() {
}

static void G__cpp_setup_func18() {
}

static void G__cpp_setup_func19() {
}

static void G__cpp_setup_func20() {
}

static void G__cpp_setup_func21() {
}

static void G__cpp_setup_func22() {
}

static void G__cpp_setup_func23() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcBosDict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
  G__cpp_setup_func13();
  G__cpp_setup_func14();
  G__cpp_setup_func15();
  G__cpp_setup_func16();
  G__cpp_setup_func17();
  G__cpp_setup_func18();
  G__cpp_setup_func19();
  G__cpp_setup_func20();
  G__cpp_setup_func21();
  G__cpp_setup_func22();
  G__cpp_setup_func23();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__BosDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__BosDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__BosDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__BosDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__BosDictLN_vectorlEunsignedsPintcOallocatorlEunsignedsPintgRsPgR = { "vector<unsigned int,allocator<unsigned int> >" , 99 , -1 };
G__linked_taginfo G__BosDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__BosDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__BosDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__BosDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__BosDictLN_TRootBeer = { "TRootBeer" , 99 , -1 };
G__linked_taginfo G__BosDictLN_RecHead_t = { "RecHead_t" , 115 , -1 };
G__linked_taginfo G__BosDictLN_recSegHead_t = { "recSegHead_t" , 115 , -1 };
G__linked_taginfo G__BosDictLN_dataSegHead_t = { "dataSegHead_t" , 115 , -1 };
G__linked_taginfo G__BosDictLN_eventInf_t = { "eventInf_t" , 115 , -1 };
G__linked_taginfo G__BosDictLN_TThread = { "TThread" , 99 , -1 };
G__linked_taginfo G__BosDictLN_vectorlETStringcOallocatorlETStringgRsPgR = { "vector<TString,allocator<TString> >" , 99 , -1 };
G__linked_taginfo G__BosDictLN_reverse_iteratorlEvectorlETStringcOallocatorlETStringgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TString,allocator<TString> >::iterator>" , 99 , -1 };
G__linked_taginfo G__BosDictLN_TBos = { "TBos" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableBosDict() {
  G__BosDictLN_TClass.tagnum = -1 ;
  G__BosDictLN_TBuffer.tagnum = -1 ;
  G__BosDictLN_TMemberInspector.tagnum = -1 ;
  G__BosDictLN_TObject.tagnum = -1 ;
  G__BosDictLN_vectorlEunsignedsPintcOallocatorlEunsignedsPintgRsPgR.tagnum = -1 ;
  G__BosDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__BosDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__BosDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__BosDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__BosDictLN_TRootBeer.tagnum = -1 ;
  G__BosDictLN_RecHead_t.tagnum = -1 ;
  G__BosDictLN_recSegHead_t.tagnum = -1 ;
  G__BosDictLN_dataSegHead_t.tagnum = -1 ;
  G__BosDictLN_eventInf_t.tagnum = -1 ;
  G__BosDictLN_TThread.tagnum = -1 ;
  G__BosDictLN_vectorlETStringcOallocatorlETStringgRsPgR.tagnum = -1 ;
  G__BosDictLN_reverse_iteratorlEvectorlETStringcOallocatorlETStringgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__BosDictLN_TBos.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableBosDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__BosDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__BosDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__BosDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__BosDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__BosDictLN_vectorlEunsignedsPintcOallocatorlEunsignedsPintgRsPgR);
   G__get_linked_tagnum_fwd(&G__BosDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__BosDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__BosDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__BosDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__BosDictLN_TRootBeer);
   G__get_linked_tagnum_fwd(&G__BosDictLN_RecHead_t);
   G__get_linked_tagnum_fwd(&G__BosDictLN_recSegHead_t);
   G__get_linked_tagnum_fwd(&G__BosDictLN_dataSegHead_t);
   G__get_linked_tagnum_fwd(&G__BosDictLN_eventInf_t);
   G__get_linked_tagnum_fwd(&G__BosDictLN_TThread);
   G__get_linked_tagnum_fwd(&G__BosDictLN_vectorlETStringcOallocatorlETStringgRsPgR);
   G__get_linked_tagnum_fwd(&G__BosDictLN_reverse_iteratorlEvectorlETStringcOallocatorlETStringgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__BosDictLN_TBos),sizeof(TBos),-1,62464,"Bos structure",G__setup_memvarTBos,G__setup_memfuncTBos);
}
extern "C" void G__cpp_setupBosDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupBosDict()");
  G__set_cpp_environmentBosDict();
  G__cpp_setup_tagtableBosDict();

  G__cpp_setup_inheritanceBosDict();

  G__cpp_setup_typetableBosDict();

  G__cpp_setup_memvarBosDict();

  G__cpp_setup_memfuncBosDict();
  G__cpp_setup_globalBosDict();
  G__cpp_setup_funcBosDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncBosDict();
  return;
}
class G__cpp_setup_initBosDict {
  public:
    G__cpp_setup_initBosDict() { G__add_setup_func("BosDict",(G__incsetup)(&G__cpp_setupBosDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initBosDict() { G__remove_setup_func("BosDict"); }
};
G__cpp_setup_initBosDict G__cpp_setup_initializerBosDict;

