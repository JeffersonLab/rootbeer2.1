//
// File generated by /apps/root/5.34.36/root//bin/rootcint at Mon Nov 25 15:22:50 2019

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME srcdIDSTWriterDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "DSTWriterDict.h"

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
   void TDSTWriter_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void delete_TDSTWriter(void *p);
   static void deleteArray_TDSTWriter(void *p);
   static void destruct_TDSTWriter(void *p);
   static void streamer_TDSTWriter(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static ROOT::TGenericClassInfo *GenerateInitInstanceLocal(const ::TDSTWriter*)
   {
      ::TDSTWriter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TDSTWriter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TDSTWriter", ::TDSTWriter::Class_Version(), "./include/TDSTWriter.h", 15,
                  typeid(::TDSTWriter), ::ROOT::DefineBehavior(ptr, ptr),
                  &::TDSTWriter::Dictionary, isa_proxy, 0,
                  sizeof(::TDSTWriter) );
      instance.SetDelete(&delete_TDSTWriter);
      instance.SetDeleteArray(&deleteArray_TDSTWriter);
      instance.SetDestructor(&destruct_TDSTWriter);
      instance.SetStreamerFunc(&streamer_TDSTWriter);
      return &instance;
   }
   ROOT::TGenericClassInfo *GenerateInitInstance(const ::TDSTWriter*)
   {
      return GenerateInitInstanceLocal((::TDSTWriter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TDSTWriter*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOTDict

//______________________________________________________________________________
atomic_TClass_ptr TDSTWriter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TDSTWriter::Class_Name()
{
   return "TDSTWriter";
}

//______________________________________________________________________________
const char *TDSTWriter::ImplFileName()
{
   return ::ROOTDict::GenerateInitInstanceLocal((const ::TDSTWriter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TDSTWriter::ImplFileLine()
{
   return ::ROOTDict::GenerateInitInstanceLocal((const ::TDSTWriter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void TDSTWriter::Dictionary()
{
   fgIsA = ::ROOTDict::GenerateInitInstanceLocal((const ::TDSTWriter*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *TDSTWriter::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gCINTMutex); if(!fgIsA) {fgIsA = ::ROOTDict::GenerateInitInstanceLocal((const ::TDSTWriter*)0x0)->GetClass();} }
   return fgIsA;
}

//______________________________________________________________________________
void TDSTWriter::Streamer(TBuffer &R__b)
{
   // Stream an object of class TDSTWriter.

   TObject::Streamer(R__b);
}

//______________________________________________________________________________
void TDSTWriter::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class TDSTWriter.
      TClass *R__cl = ::TDSTWriter::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fRootbeer", &fRootbeer);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*bankTree", &bankTree);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fRootFile", &fRootFile);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fBrName[100]", fBrName);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fBrFormat[100]", fBrFormat);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fBigBuff", &fBigBuff);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fBuffPtr", &fBuffPtr);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fUsedBankTot", &fUsedBankTot);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fUsedBankIndex[220]", fUsedBankIndex);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fUsedBankRows[220]", fUsedBankRows);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fUsedBankType[220]", fUsedBankType);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fUsedBankName[220][5]", fUsedBankName);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fDropBankName[220][5]", fDropBankName);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fDropBankTotal", &fDropBankTotal);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fDropBankFlag", &fDropBankFlag);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fDSTInitFlag", &fDSTInitFlag);
      TObject::ShowMembers(R__insp);
}

namespace ROOTDict {
   // Wrapper around operator delete
   static void delete_TDSTWriter(void *p) {
      delete ((::TDSTWriter*)p);
   }
   static void deleteArray_TDSTWriter(void *p) {
      delete [] ((::TDSTWriter*)p);
   }
   static void destruct_TDSTWriter(void *p) {
      typedef ::TDSTWriter current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TDSTWriter(TBuffer &buf, void *obj) {
      ((::TDSTWriter*)obj)->::TDSTWriter::Streamer(buf);
   }
} // end of namespace ROOTDict for class ::TDSTWriter

/********************************************************
* src/DSTWriterDict.cxx
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

extern "C" void G__cpp_reset_tagtableDSTWriterDict();

extern "C" void G__set_cpp_environmentDSTWriterDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("include/TDSTWriter.h");
  G__cpp_reset_tagtableDSTWriterDict();
}
#include <new>
extern "C" int G__cpp_dllrevDSTWriterDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* TDSTWriter */
static int G__DSTWriterDict_775_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TDSTWriter* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new TDSTWriter((TFile*) G__int(libp->para[0]));
   } else {
     p = new((void*) gvp) TDSTWriter((TFile*) G__int(libp->para[0]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__DSTWriterDictLN_TDSTWriter));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DSTWriterDict_775_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((TDSTWriter*) G__getstructoffset())->Init((TRootBeer*) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DSTWriterDict_775_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((TDSTWriter*) G__getstructoffset())->WriteDST());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DSTWriterDict_775_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((TDSTWriter*) G__getstructoffset())->DropBank((char*) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DSTWriterDict_775_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((const TDSTWriter*) G__getstructoffset())->GetTree());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DSTWriterDict_775_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) TDSTWriter::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DSTWriterDict_775_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TDSTWriter::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DSTWriterDict_775_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) TDSTWriter::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DSTWriterDict_775_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TDSTWriter::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DSTWriterDict_775_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TDSTWriter*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DSTWriterDict_775_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TDSTWriter::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DSTWriterDict_775_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TDSTWriter::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DSTWriterDict_775_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TDSTWriter::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DSTWriterDict_775_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TDSTWriter::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__DSTWriterDict_775_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   TDSTWriter* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new TDSTWriter(*(TDSTWriter*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__DSTWriterDictLN_TDSTWriter));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef TDSTWriter G__TTDSTWriter;
static int G__DSTWriterDict_775_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (TDSTWriter*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((TDSTWriter*) (soff+(sizeof(TDSTWriter)*i)))->~G__TTDSTWriter();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (TDSTWriter*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((TDSTWriter*) (soff))->~G__TTDSTWriter();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__DSTWriterDict_775_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TDSTWriter* dest = (TDSTWriter*) G__getstructoffset();
   *dest = *(TDSTWriter*) libp->para[0].ref;
   const TDSTWriter& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* TDSTWriter */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncDSTWriterDict {
 public:
  G__Sizep2memfuncDSTWriterDict(): p(&G__Sizep2memfuncDSTWriterDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncDSTWriterDict::*p)();
};

size_t G__get_sizep2memfuncDSTWriterDict()
{
  G__Sizep2memfuncDSTWriterDict a;
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
extern "C" void G__cpp_setup_inheritanceDSTWriterDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__DSTWriterDictLN_TDSTWriter))) {
     TDSTWriter *G__Lderived;
     G__Lderived=(TDSTWriter*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__DSTWriterDictLN_TDSTWriter),G__get_linked_tagnum(&G__DSTWriterDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableDSTWriterDict() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DSTWriterDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DSTWriterDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DSTWriterDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DSTWriterDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<UInt_t>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_vectorlEunsignedsPintcOallocatorlEunsignedsPintgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TString>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_vectorlETStringcOallocatorlETStringgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_reverse_iteratorlEvectorlETStringcOallocatorlETStringgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DSTWriterDictLN_vectorlETStringcOallocatorlETStringgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_reverse_iteratorlEvectorlETStringcOallocatorlETStringgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DSTWriterDictLN_vectorlETStringcOallocatorlETStringgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<std::bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long,const TObject**>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<std::string,TObjArray*>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<string,TObjArray*>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<string,TObjArray*>",117,G__get_linked_tagnum(&G__DSTWriterDictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<string,TObjArray*,less<string> >",117,G__get_linked_tagnum(&G__DSTWriterDictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* TDSTWriter */
static void G__setup_memvarTDSTWriter(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__DSTWriterDictLN_TDSTWriter));
   { TDSTWriter *p; p=(TDSTWriter*)0x1000; if (p) { }
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__DSTWriterDictLN_TRootBeer),-1,-1,2,"fRootbeer=",0,"the main object");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__DSTWriterDictLN_TTree),-1,-1,2,"bankTree=",0,"Tree to hold banks ");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__DSTWriterDictLN_TFile),-1,-1,2,"fRootFile=",0,"root file device");
   G__memvar_setup((void*)0,99,0,0,-1,-1,-1,2,"fBrName[100]=",0,"branch name");
   G__memvar_setup((void*)0,99,0,0,-1,-1,-1,2,"fBrFormat[100]=",0,"branch format string");
   G__memvar_setup((void*)0,67,0,0,-1,-1,-1,2,"fBigBuff=",0,"allocate a 100k buffer for multi instance banks");
   G__memvar_setup((void*)0,67,0,0,-1,-1,-1,2,"fBuffPtr=",0,"running pointer within bigBuff");
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,2,"fUsedBankTot=",0,"total no of banks accessed");
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,2,"fUsedBankIndex[220]=",0,"indexes in bankAddress[]			");
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,2,"fUsedBankRows[220]=",0,"total no of rows for mulit instance banks");
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,2,"fUsedBankType[220]=",0,"0 for singles 1 for multis");
   G__memvar_setup((void*)0,99,0,0,-1,-1,-1,2,"fUsedBankName[220][5]=",0,"names of used banks");
   G__memvar_setup((void*)0,99,0,0,-1,-1,-1,2,"fDropBankName[220][5]=",0,"names of banks to be dropped (ie not written to DST )");
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,2,"fDropBankTotal=",0,"counter for no od dropped banks");
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,2,"fDropBankFlag=",0,(char*)NULL);
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,2,"fDSTInitFlag=",0,"Flag if dst initislised");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__DSTWriterDictLN_TClass),G__defined_typename("atomic_TClass_ptr"),-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarDSTWriterDict() {
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
static void G__setup_memfuncTDSTWriter(void) {
   /* TDSTWriter */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__DSTWriterDictLN_TDSTWriter));
   G__memfunc_setup("TDSTWriter",956,G__DSTWriterDict_775_0_1, 105, G__get_linked_tagnum(&G__DSTWriterDictLN_TDSTWriter), -1, 0, 1, 1, 1, 0, "U 'TFile' - 0 - -", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Init",404,G__DSTWriterDict_775_0_2, 105, -1, -1, 0, 1, 1, 1, 0, "U 'TRootBeer' - 0 - -", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("WriteDST",758,G__DSTWriterDict_775_0_3, 105, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DropBank",785,G__DSTWriterDict_775_0_4, 105, -1, -1, 0, 1, 1, 1, 0, "C - - 0 - -", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetTree",688,G__DSTWriterDict_775_0_5, 85, G__get_linked_tagnum(&G__DSTWriterDictLN_TTree), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__DSTWriterDict_775_0_6, 85, G__get_linked_tagnum(&G__DSTWriterDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&TDSTWriter::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__DSTWriterDict_775_0_7, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TDSTWriter::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__DSTWriterDict_775_0_8, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&TDSTWriter::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__DSTWriterDict_775_0_9, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&TDSTWriter::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__DSTWriterDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__DSTWriterDict_775_0_13, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__DSTWriterDict_775_0_14, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TDSTWriter::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__DSTWriterDict_775_0_15, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TDSTWriter::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__DSTWriterDict_775_0_16, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TDSTWriter::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__DSTWriterDict_775_0_17, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TDSTWriter::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("TDSTWriter", 956, G__DSTWriterDict_775_0_18, (int) ('i'), G__get_linked_tagnum(&G__DSTWriterDictLN_TDSTWriter), -1, 0, 1, 1, 1, 0, "u 'TDSTWriter' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~TDSTWriter", 1082, G__DSTWriterDict_775_0_19, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__DSTWriterDict_775_0_20, (int) ('u'), G__get_linked_tagnum(&G__DSTWriterDictLN_TDSTWriter), -1, 1, 1, 1, 1, 0, "u 'TDSTWriter' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncDSTWriterDict() {
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
}

static void G__cpp_setup_global4() {
}

static void G__cpp_setup_global5() {
}

static void G__cpp_setup_global6() {
}

static void G__cpp_setup_global7() {
}

static void G__cpp_setup_global8() {
}

static void G__cpp_setup_global9() {
}

static void G__cpp_setup_global10() {
}

static void G__cpp_setup_global11() {
}

static void G__cpp_setup_global12() {
}

static void G__cpp_setup_global13() {
}

static void G__cpp_setup_global14() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalDSTWriterDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
  G__cpp_setup_global2();
  G__cpp_setup_global3();
  G__cpp_setup_global4();
  G__cpp_setup_global5();
  G__cpp_setup_global6();
  G__cpp_setup_global7();
  G__cpp_setup_global8();
  G__cpp_setup_global9();
  G__cpp_setup_global10();
  G__cpp_setup_global11();
  G__cpp_setup_global12();
  G__cpp_setup_global13();
  G__cpp_setup_global14();
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
}

static void G__cpp_setup_func24() {
}

static void G__cpp_setup_func25() {
}

static void G__cpp_setup_func26() {
}

static void G__cpp_setup_func27() {
}

static void G__cpp_setup_func28() {
}

static void G__cpp_setup_func29() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcDSTWriterDict() {
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
  G__cpp_setup_func24();
  G__cpp_setup_func25();
  G__cpp_setup_func26();
  G__cpp_setup_func27();
  G__cpp_setup_func28();
  G__cpp_setup_func29();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__DSTWriterDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__DSTWriterDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__DSTWriterDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__DSTWriterDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__DSTWriterDictLN_vectorlEunsignedsPintcOallocatorlEunsignedsPintgRsPgR = { "vector<unsigned int,allocator<unsigned int> >" , 99 , -1 };
G__linked_taginfo G__DSTWriterDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__DSTWriterDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__DSTWriterDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__DSTWriterDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__DSTWriterDictLN_TRootBeer = { "TRootBeer" , 99 , -1 };
G__linked_taginfo G__DSTWriterDictLN_vectorlETStringcOallocatorlETStringgRsPgR = { "vector<TString,allocator<TString> >" , 99 , -1 };
G__linked_taginfo G__DSTWriterDictLN_reverse_iteratorlEvectorlETStringcOallocatorlETStringgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TString,allocator<TString> >::iterator>" , 99 , -1 };
G__linked_taginfo G__DSTWriterDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR = { "iterator<bidirectional_iterator_tag,TObject*,long,const TObject**,const TObject*&>" , 115 , -1 };
G__linked_taginfo G__DSTWriterDictLN_TTree = { "TTree" , 99 , -1 };
G__linked_taginfo G__DSTWriterDictLN_TFile = { "TFile" , 99 , -1 };
G__linked_taginfo G__DSTWriterDictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR = { "map<string,TObjArray*,less<string>,allocator<pair<const string,TObjArray*> > >" , 99 , -1 };
G__linked_taginfo G__DSTWriterDictLN_TDSTWriter = { "TDSTWriter" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableDSTWriterDict() {
  G__DSTWriterDictLN_TClass.tagnum = -1 ;
  G__DSTWriterDictLN_TBuffer.tagnum = -1 ;
  G__DSTWriterDictLN_TMemberInspector.tagnum = -1 ;
  G__DSTWriterDictLN_TObject.tagnum = -1 ;
  G__DSTWriterDictLN_vectorlEunsignedsPintcOallocatorlEunsignedsPintgRsPgR.tagnum = -1 ;
  G__DSTWriterDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__DSTWriterDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__DSTWriterDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__DSTWriterDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__DSTWriterDictLN_TRootBeer.tagnum = -1 ;
  G__DSTWriterDictLN_vectorlETStringcOallocatorlETStringgRsPgR.tagnum = -1 ;
  G__DSTWriterDictLN_reverse_iteratorlEvectorlETStringcOallocatorlETStringgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__DSTWriterDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR.tagnum = -1 ;
  G__DSTWriterDictLN_TTree.tagnum = -1 ;
  G__DSTWriterDictLN_TFile.tagnum = -1 ;
  G__DSTWriterDictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR.tagnum = -1 ;
  G__DSTWriterDictLN_TDSTWriter.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableDSTWriterDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_vectorlEunsignedsPintcOallocatorlEunsignedsPintgRsPgR);
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_TRootBeer);
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_vectorlETStringcOallocatorlETStringgRsPgR);
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_reverse_iteratorlEvectorlETStringcOallocatorlETStringgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR);
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_TTree);
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_TFile);
   G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__DSTWriterDictLN_TDSTWriter),sizeof(TDSTWriter),-1,62464,(char*)NULL,G__setup_memvarTDSTWriter,G__setup_memfuncTDSTWriter);
}
extern "C" void G__cpp_setupDSTWriterDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupDSTWriterDict()");
  G__set_cpp_environmentDSTWriterDict();
  G__cpp_setup_tagtableDSTWriterDict();

  G__cpp_setup_inheritanceDSTWriterDict();

  G__cpp_setup_typetableDSTWriterDict();

  G__cpp_setup_memvarDSTWriterDict();

  G__cpp_setup_memfuncDSTWriterDict();
  G__cpp_setup_globalDSTWriterDict();
  G__cpp_setup_funcDSTWriterDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncDSTWriterDict();
  return;
}
class G__cpp_setup_initDSTWriterDict {
  public:
    G__cpp_setup_initDSTWriterDict() { G__add_setup_func("DSTWriterDict",(G__incsetup)(&G__cpp_setupDSTWriterDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initDSTWriterDict() { G__remove_setup_func("DSTWriterDict"); }
};
G__cpp_setup_initDSTWriterDict G__cpp_setup_initializerDSTWriterDict;

