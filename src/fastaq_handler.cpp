#include <string>
#include <cstring>
#include <iostream>
//#include <fstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "fastaq_handler.h"
#include "utils.h"

FastaqHandler::FastaqHandler(const std::string& filepath)
    : closed_status(3)
    , filepath(filepath)
    , gzipped(false)
    , num_reads_parsed(0)
    , read_status(0)
{
    // level for boost logging
    //    logging::core::get()->set_filter(logging::trivial::severity >= g_log_level);

    BOOST_LOG_TRIVIAL(debug) << "Open fastaq file" << filepath;
    this->fastaq_file = gzopen(filepath.c_str(), "r");
    if (this->fastaq_file == nullptr) {
        throw "Unable to open " + this->filepath;
    }
    this->inbuf = kseq_init(this->fastaq_file);
}

FastaqHandler::~FastaqHandler() {
    if (!this->is_closed()) {
        close();
    }
}

bool FastaqHandler::eof() const { return (this->read_status == -1); }

void FastaqHandler::get_next()
{
    this->read_status = kseq_read(this->inbuf);

    if (this->eof()) {
        return;
    }

    if (this->read_status == -2) {
        throw "Truncated quality string detected";
    } else if (this->read_status == -3) {
        throw "Error reading " + this->filepath;
    }

    ++this->num_reads_parsed;
    this->name = this->inbuf->name.s;
    this->read = this->inbuf->seq.s;
}

void FastaqHandler::skip_next()
{
    if (this->eof()) {
        return;
    }
    this->get_next();
}

void print(std::ifstream& infile)
{
    char file;
    std::vector<char> read;
    uint i = 0;

    // Read infile to vector
    while (!infile.eof()) {
        infile >> file;
        read.push_back(file);
    }

    // Print read vector
    for (i = 0; i < read.size(); i++) {
        std::cout << read[i];
    }
}

void print(std::istream& infile)
{
    char file;
    std::vector<char> read;
    int i = 0;

    // Read infile to vector
    while (!infile.eof()) {
        infile >> file;
        read.push_back(file);
    }

    // Print read vector
    for (auto i = 0; i < read.size(); i++) {
        std::cout << read[i];
    }
}

void FastaqHandler::get_id(const uint32_t& id)
{
    const uint32_t one_based_id = id + 1;
    if (one_based_id < this->num_reads_parsed) {
        BOOST_LOG_TRIVIAL(warning)
            << "restart buffer as have id " << num_reads_parsed << " and want id "
            << one_based_id << " (" << id << ") with 0-based indexing.";
        num_reads_parsed = 0;
        name.clear();
        read.clear();
        gzrewind(this->fastaq_file);
        kseq_rewind(this->inbuf);
    }

    while (id > 1 and num_reads_parsed < id) {
        skip_next();
        if (eof()) {
            break;
        }
    }

    while (num_reads_parsed <= id) {
        get_next();
        if (eof() && num_reads_parsed < id) {
            throw std::out_of_range("Requested a read past the end of the file.");
        }
    }
}

void FastaqHandler::close()
{
    BOOST_LOG_TRIVIAL(debug) << "Close fastaq file";
    this->closed_status = gzclose(this->fastaq_file);
}

bool FastaqHandler::is_closed() const { return this->closed_status == Z_OK; }