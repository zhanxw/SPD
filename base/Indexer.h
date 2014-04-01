#ifndef _INDEXER_H_
#define _INDEXER_H_

#include <map>
#include <set>
#include <string>
#include <vector>

/**
 * Indexer class is convenient tool to index a vector class
 */
class Indexer{
 public:
  Indexer(const std::vector<std::string>& a) {
    for (size_t i = 0; i < a.size(); ++i) {
      if (m.count(a[i]) == 0) {
        m[a[i]] = i;
      } else {
        duplication.insert(a[i]);
      }
    }
  }
  bool hasDuplication() {
    return this->duplication.size() > 0;
  }
  int operator[](const std::string& s) const {
    if (m.count(s) == 0) return -1;
    return m.find(s)->second;
  }
 private:
  std::map<std::string, int> m;
  std::set<std::string> duplication;
};


#endif /* _INDEXER_H_ */
