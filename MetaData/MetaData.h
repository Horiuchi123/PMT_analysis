#ifndef RootInterface_MetaData_h
#define RootInterface_MetaData_h

#include "TObject.h"

#include <string>
#include <any>
#include <map>

class MetaData : public TObject {

public:
  /**
   * @brief Method that prints current metadata entries */
  void Print();

  /**
   * @brief Method that checks if this instance has a key */
  bool HasKey(std::string key) { return Map.count(key)>0 ? true : false; }

  /**
   * @brief For a bit more easy access to the data */
  std::string& operator[](const std::string key) { return Map[key]; }

  /**
   * @brief Actual data container */
  std::map<std::string,std::string> Map;

  ClassDef(MetaData,1);
};

#endif
