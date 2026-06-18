#include "syncfile.hpp"
#include <iostream>

SyncFile::SyncFile(const std::string& path) {
  std::cout
      << "Constructing SyncFile at "
      << this
      << std::endl;

  _csv_output.open(path);
}

SyncFile::~SyncFile() {

  std::cout
      << "Rows written = "
      << rows_written.load()
      << std::endl;

  _csv_output.close();
}

bool SyncFile::write(const std::string& data) {
  try {
    std::lock_guard<std::mutex> lock(_writerMutex);
    _csv_output << data << "\n";
    return true;
  } catch (...) {
    return false;
  }
}

bool SyncFile::write(const csv_data& data) {
  try {
    std::lock_guard<std::mutex> lock(_writerMutex);

    _csv_output << data << "\n";

    rows_written.fetch_add(1);

    if (!_csv_output.good()) {
      std::cerr << "WRITE FAILED" << std::endl;
      return false;
    }

    return true;

  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return false;
  }
}

// #include "syncfile.hpp"

// // SyncFile::SyncFile(const std::string& path) : _path(path) { _csv_output.open(_path); }
// SyncFile::SyncFile(const std::string& path) {
//   _csv_output.open(path);
// }
// // SyncFile::~SyncFile() {
// //   writeToFile();
// //   _csv_output.close();
// // }
// SyncFile::~SyncFile() {
//   _csv_output.close();
// }

// bool SyncFile::writeToFile() {
//   try {
//     // Get the file lock for this thread
//     std::lock_guard<std::mutex> lock(_writerMutex);
//     // loop over queue wring out events to the file before closing
//     // while (!_writeQueue.empty()) {
//     //   _csv_output << _writeQueue.front() << "\n";
//     //   _writeQueue.pop();
//     // }
//     return true;
//   } catch (const std::exception& e) {
//     std::cerr << e.what() << '\n';
//     return false;
//   }
// }

// bool SyncFile::write(const std::string& data) {
//   try {
//     std::lock_guard<std::mutex> lock(_writerMutex);
//     _csv_output << data << std::endl;
//     return true;
//   } catch (const std::exception& e) {
//     std::cerr << e.what() << '\n';
//     return false;
//   }
// }

// // bool SyncFile::write(const csv_data& data) {
// //   try {
// //     // Get the file lock for this thread
// //     std::lock_guard<std::mutex> lock(_writerMutex);
// //     // put data into a queue to be written later
// //     // _writeQueue.push(data);
// //     //_csv_output << data << std::endl;
// //     return true;
// //   } catch (const std::exception& e) {
// //     std::cerr << e.what() << '\n';
// //     return false;
// //   }
// // }

// bool SyncFile::write(const csv_data& data) {
//   try {
//     std::lock_guard<std::mutex> lock(_writerMutex);

//     _csv_output
//       << data.run << ","
//       << data.event << ","
//       << data.pid_prot_rec << ","
//       << data.pid_pip_rec << ","
//       << data.pid_pim_rec << ","
//       << data.w << ","
//       << data.q2 << "\n";

//     return true;
//   } catch (...) {
//     return false;
//   }
// }

// // bool SyncFile::write(const csv_data& data) {
// //   try {
// //     std::lock_guard<std::mutex> lock(_writerMutex);
// //     _csv_output << data << '\n';
// //     return true;
// //   } catch (const std::exception& e) {
// //     std::cerr << e.what() << '\n';
// //     return false;
// //   }
// // }
