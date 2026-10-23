// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME MetaDataDict
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
#include "MetaData.h"
#include "MetaData.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_MetaData(void *p = nullptr);
   static void *newArray_MetaData(Long_t size, void *p);
   static void delete_MetaData(void *p);
   static void deleteArray_MetaData(void *p);
   static void destruct_MetaData(void *p);
   static void streamer_MetaData(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MetaData*)
   {
      ::MetaData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MetaData >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("MetaData", ::MetaData::Class_Version(), "MetaData.h", 10,
                  typeid(::MetaData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MetaData::Dictionary, isa_proxy, 16,
                  sizeof(::MetaData) );
      instance.SetNew(&new_MetaData);
      instance.SetNewArray(&newArray_MetaData);
      instance.SetDelete(&delete_MetaData);
      instance.SetDeleteArray(&deleteArray_MetaData);
      instance.SetDestructor(&destruct_MetaData);
      instance.SetStreamerFunc(&streamer_MetaData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MetaData*)
   {
      return GenerateInitInstanceLocal(static_cast<::MetaData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::MetaData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr MetaData::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *MetaData::Class_Name()
{
   return "MetaData";
}

//______________________________________________________________________________
const char *MetaData::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MetaData*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int MetaData::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MetaData*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MetaData::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MetaData*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MetaData::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MetaData*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void MetaData::Streamer(TBuffer &R__b)
{
   // Stream an object of class MetaData.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      {
         map<string,string> &R__stl =  Map;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         for (R__i = 0; R__i < R__n; R__i++) {
            string R__t;
            {TString R__str;
             R__str.Streamer(R__b);
             R__t = R__str.Data();}
            string R__t2;
            {TString R__str;
             R__str.Streamer(R__b);
             R__t2 = R__str.Data();}
            typedef class std::basic_string<char> Value_t;
            std::pair<Value_t const, class std::basic_string<char> > R__t3(R__t,R__t2);
            R__stl.insert(R__t3);
         }
      }
      R__b.CheckByteCount(R__s, R__c, MetaData::IsA());
   } else {
      R__c = R__b.WriteVersion(MetaData::IsA(), kTRUE);
      TObject::Streamer(R__b);
      {
         map<string,string> &R__stl =  Map;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            map<string,string>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            {TString R__str(((*R__k).first ).c_str());
             R__str.Streamer(R__b);};
            {TString R__str(((*R__k).second).c_str());
             R__str.Streamer(R__b);};
            }
         }
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MetaData(void *p) {
      return  p ? new(p) ::MetaData : new ::MetaData;
   }
   static void *newArray_MetaData(Long_t nElements, void *p) {
      return p ? new(p) ::MetaData[nElements] : new ::MetaData[nElements];
   }
   // Wrapper around operator delete
   static void delete_MetaData(void *p) {
      delete (static_cast<::MetaData*>(p));
   }
   static void deleteArray_MetaData(void *p) {
      delete [] (static_cast<::MetaData*>(p));
   }
   static void destruct_MetaData(void *p) {
      typedef ::MetaData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_MetaData(TBuffer &buf, void *obj) {
      ((::MetaData*)obj)->::MetaData::Streamer(buf);
   }
} // end of namespace ROOT for class ::MetaData

namespace ROOT {
   static TClass *maplEstringcOstringgR_Dictionary();
   static void maplEstringcOstringgR_TClassManip(TClass*);
   static void *new_maplEstringcOstringgR(void *p = nullptr);
   static void *newArray_maplEstringcOstringgR(Long_t size, void *p);
   static void delete_maplEstringcOstringgR(void *p);
   static void deleteArray_maplEstringcOstringgR(void *p);
   static void destruct_maplEstringcOstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,string>*)
   {
      map<string,string> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,string>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,string>", -2, "map", 102,
                  typeid(map<string,string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOstringgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,string>) );
      instance.SetNew(&new_maplEstringcOstringgR);
      instance.SetNewArray(&newArray_maplEstringcOstringgR);
      instance.SetDelete(&delete_maplEstringcOstringgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOstringgR);
      instance.SetDestructor(&destruct_maplEstringcOstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,string> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,string>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,string>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,string>*>(nullptr))->GetClass();
      maplEstringcOstringgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOstringgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,string> : new map<string,string>;
   }
   static void *newArray_maplEstringcOstringgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,string>[nElements] : new map<string,string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOstringgR(void *p) {
      delete (static_cast<map<string,string>*>(p));
   }
   static void deleteArray_maplEstringcOstringgR(void *p) {
      delete [] (static_cast<map<string,string>*>(p));
   }
   static void destruct_maplEstringcOstringgR(void *p) {
      typedef map<string,string> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,string>

namespace ROOT {
   // Registration Schema evolution read functions
   int RecordReadRules_MetaDataDict() {
      return 0;
   }
   static int _R__UNIQUE_DICT_(ReadRules_MetaDataDict) = RecordReadRules_MetaDataDict();R__UseDummy(_R__UNIQUE_DICT_(ReadRules_MetaDataDict));
} // namespace ROOT
namespace {
  void TriggerDictionaryInitialization_MetaDataDict_Impl() {
    static const char* headers[] = {
"MetaData.h",
nullptr
    };
    static const char* includePaths[] = {
"/home/hkpd/root/include/",
"/home/hkpd/individuals/horiuchi/robo/analysis/cpp/MetaData/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MetaDataDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$MetaData.h")))  MetaData;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MetaDataDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "MetaData.h"
//#include "RootInterface/TH1Util.h"
//#include "RootInterface/MetaData.h"
#include "MetaData.h"
//#include "RootInterface/Hit.h"
//#include <any>
#ifdef __CINT__
//#pragma link off all globals;
//#pragma link off all classes;
//#pragma link off all functions;
//#pragma link C++ nestedclass;
//#pragma link C++ namespace TH1Util;
//#pragma link C++ function TH1Util::ReduceHistAxis;
//#pragma link C++ function TH1Util::IntegratedError;
//#pragma link C++ function TH1Util::GetLogarithmicYmin;
//#pragma link C++ function TH1Util::GetLogarithmicYmax;
//#pragma link C++ function TH1Util::SetHistTitleSize;
//#pragma link C++ function TH1Util::HistScaleX;
//#pragma link C++ function SetAtlasStyle;
#pragma link C++ class MetaData;
//#pragma link C++ class Hit;
//#pragma link C++ class std::vector<Hit>;
#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"MetaData", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MetaDataDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MetaDataDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MetaDataDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MetaDataDict() {
  TriggerDictionaryInitialization_MetaDataDict_Impl();
}
