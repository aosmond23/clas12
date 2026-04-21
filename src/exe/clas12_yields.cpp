#include "clas12_yields.hpp"
#include <future>
#include <thread>

// Define the static member
bool csv_data::isGenerated = true;  // Default to true (Generated data)

int main(int argc, char** argv) {
        // Need this to make sure root doesn't break
        ROOT::EnableThreadSafety();
        // std::ios::sync_with_stdio(false);

        // Make sure we don't create more threads than files
        int NUM_THREADS = 4;
        if (getenv("NUM_THREADS") != NULL) NUM_THREADS = atoi(getenv("NUM_THREADS"));
        // if (NUM_THREADS > argc - NUM_THREADS) NUM_THREADS = 1;

        int num_inputs = argc - 2; // argv[0]=program, argv[1]=output file

        if (NUM_THREADS > num_inputs) {NUM_THREADS = num_inputs;}

        if (NUM_THREADS < 1) {NUM_THREADS = 1;}

        // Make a vector of vectors of strings the size of the number of threads
        std::vector<std::vector<std::string> > infilenames(NUM_THREADS);
        // Get the output file name
        std::string outfilename;

        if (argc >= 2) {
                // First argument is the output file
                outfilename = argv[1];
                // Set the output type based on the output filename; 9/5/24
                csv_data::setOutputType(outfilename);
                // All other files are split evently by the under of threads
                // for (int i = 2; i < argc; i++) infilenames[i % NUM_THREADS].push_back(argv[i]);
                int t = 0;
                for (int i = 2; i < argc; i++) {infilenames[t].push_back(argv[i]); t = (t + 1) % NUM_THREADS;}
        } else {
                return 1;
        }

        // Make your histograms object as a shared pointer that all the threads will have
        auto csv_output_file = std::make_shared<SyncFile>(outfilename);
        csv_output_file->write(csv_data::header());
        // auto run_files = [&csv_output_file](auto&& inputs, auto&& thread_id) mutable {
        auto run_files = [&csv_output_file, &outfilename](std::vector<std::string> inputs, auto&& thread_id) mutable {
                 // Called once for each thread
                 // Make a new chain to process for this thread
                 auto chain = std::make_shared<TChain>("clas12");

                 // Add every file to the chain
                 for (auto in : inputs) chain->Add(in.c_str());
                                
                 // Run the function over each thread
                 // return run(chain, csv_output_file, thread_id);

                 return run<Pass2_Cuts>(std::move(chain), csv_output_file, thread_id, outfilename); // commented out 9/4/24


         };

        // Make a set of threads (Futures are special threads which return a value)
        std::future<size_t> threads[NUM_THREADS];

        // Define events to be used to get Hz later
        size_t events = 0;

        // Start timer
        auto start = std::chrono::high_resolution_clock::now();
        // For each thread
        for (size_t i = 0; i < NUM_THREADS; i++) {
                // Set the thread to run a task A-Syncroisly
                // The function we run is the first argument (run_files)
                // The functions areruments are all the remaining arguments
                threads[i] = std::async(run_files, infilenames.at(i), i);
        }

        // For each thread
        for (size_t i = 0; i < NUM_THREADS; i++) {
                // Get the information from the thread in this case how many events each thread actually computed
                events += threads[i].get();
        }

        // Timer and Hz calculator functions that print at the end
        std::cout.imbue(std::locale("")); // Puts commas in
        std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
        std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
        std::cout << BOLDYELLOW << events / elapsed_full.count() << " Hz" << DEF << std::endl;

        csv_output_file->writeToFile();

        std::chrono::duration<double> elapsed_full_write = (std::chrono::high_resolution_clock::now() - start);
        std::cout << RED << elapsed_full_write.count() << " sec" << DEF << std::endl;
        std::cout << BOLDYELLOW << events / elapsed_full_write.count() << " Hz" << DEF << std::endl;

        return 0;
}
