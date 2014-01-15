#include "MatFile.h"
#include <stdint.h>
#include "sstream"


////////////////////
// Implementation //
////////////////////

MatFile::MatFile(const std::string& fname)
{
    fs_.open(
        fname.c_str(), 
        std::ofstream::out | std::ofstream::binary | std::ofstream::trunc
    );
    writeHeader();
}

void MatFile::startCell(const std::string& fieldName, int n)
{
    int32_t header[] = {
        // Header
        14, 0,
        
        // Flags
        6, 8, 1, 0,

        // Sizes
        5, 4, n, 0,

        // Name
        1, fieldName.size()
    };

    // Name padding
    const char padding[] = "        ";
    const int npd = (fieldName.size()+7)/8*8 - fieldName.size();

    // Write to file
    fs_.write(reinterpret_cast<const char*>(header), 4);
    lengthPos_ = fs_.tellp();
    fs_.write(reinterpret_cast<const char*>(header+1), 11*4);
    fs_.write(fieldName.c_str(), fieldName.size());
    fs_.write(padding, npd);

    // Ready for writing fields...
}

void MatFile::endCell()
{
    int endPos = fs_.tellp();
    int32_t length = endPos - lengthPos_ - 4;
    fs_.seekp(lengthPos_);
    fs_.write(reinterpret_cast<const char*>(&length), 4);
    fs_.seekp(endPos);
}


void MatFile::close()
{
    fs_.close();
}

void MatFile::writeHeader()
{
    static const char description[] =
        "Stencil-parareal project. Output data file. Good luck           "
        "                                                                ";
    const uint16_t version = 0x0100;
    const uint16_t endian = 0x4D49;

    fs_.seekp(0);
    fs_.write(description, 124);
    fs_.write(reinterpret_cast<const char*>(&version), 2);
    fs_.write(reinterpret_cast<const char*>(&endian), 2);
}

