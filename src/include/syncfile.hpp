#ifndef SYNCFILE_H_GUARD
#define SYNCFILE_H_GUARD

#include <fstream>
#include <mutex>
#include "csv_data.hpp"
#include <atomic>

class SyncFile {
 public:
  SyncFile(const std::string& path);
  ~SyncFile();

  bool write(const csv_data& data);
  bool write(const std::string& data);

 private:
  std::ofstream _csv_output;
  std::mutex _writerMutex;
  std::atomic<size_t> rows_written{0};
};

#endif


// #ifndef SYNCFILE_H_GUARD
// #define SYNCFILE_H_GUARD

// #include <fstream>
// #include <iomanip>
// #include <iostream>
// #include <mutex>
// #include <queue>
// #include <string>
// #include "csv_data.hpp"

// // https://stackoverflow.com/questions/33596910/c-how-to-have-multiple-threads-write-to-a-file

// // class SyncFile {
// //  public:
// //   SyncFile(const std::string& path);
// //   ~SyncFile();
// //   bool write(const csv_data& data);
// //   bool write(const std::string& data);
// //   bool writeToFile();

// //  private:
// //   std::string _path;
// //   std::ofstream _csv_output;
// //   std::mutex _writerMutex;
// //   std::queue<csv_data> _writeQueue;
// // };

// class SyncFile {
//  public:
//   SyncFile(const std::string& path);
//   ~SyncFile();

//   bool write(const csv_data& data);

//  private:
//   std::ofstream _csv_output;
//   std::mutex _writerMutex;
// };

// #endif