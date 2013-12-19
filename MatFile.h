#include <iostream>
#include <fstream>
#include <string>
#include "SharedInfrastructure.h"

class MatFile
{
public:
    MatFile(const std::string& fname);
    void startCell(const std::string&, int n);
    void endCell();
    void addField(const IJKRealField& field, int inc = -1);

private:
    std::ofstream fs_;
    int lengthPos_;

    void writeHeader();

};
