#include "syncfile.hpp"

SyncFile::SyncFile(const std::string& path) : _path(path) { _csv_output.open(_path); }
SyncFile::~SyncFile() {
  writeToFile();
  _csv_output.close();
}

bool SyncFile::writeToFile() {
  try {
    // Get the file lock for this thread
    std::lock_guard<std::mutex> lock(_writerMutex);
    // loop over queue wring out events to the file before closing
    while (!_writeQueue.empty()) {
      _csv_output << _writeQueue.front() << "\n";
      _writeQueue.pop();
    }
    return true;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return false;
  }
}

bool SyncFile::write(const std::string& data) {
  try {
    std::lock_guard<std::mutex> lock(_writerMutex);
    _csv_output << data << std::endl;
    return true;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return false;
  }
}

// bool SyncFile::write(const csv_data& data) {
//   try {
//     // Get the file lock for this thread
//     std::lock_guard<std::mutex> lock(_writerMutex);
//     // put data into a queue to be written later
//     _writeQueue.push(data);
//     //_csv_output << data << std::endl;
//     return true;
//   } catch (const std::exception& e) {
//     std::cerr << e.what() << '\n';
//     return false;
//   }
// }

bool SyncFile::write(const csv_data& data) {
  try {
    std::lock_guard<std::mutex> lock(_writerMutex);
    _csv_output << data << '\n';
    return true;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return false;
  }
}
