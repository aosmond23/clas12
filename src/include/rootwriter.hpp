#ifndef ROOT_WRITER_H_GUARD
#define ROOT_WRITER_H_GUARD

#include <string>
#include "TFile.h"
#include "TTree.h"
#include "csv_data.hpp"

class RootWriter {

public:

    RootWriter(const std::string& filename, bool isGenerated);
    ~RootWriter();

    void Fill(const csv_data& data);

private:
    // void createGenBranches();
    // void createRecBranches();
    void createBranches();

    bool generated;

    TFile* file;
    TTree* tree;

    csv_data buffer;
};

#endif