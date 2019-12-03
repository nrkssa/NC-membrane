#ifndef _PARSER_H__
#define _PARSER_H__
#include <vector>
#include "_constant.h"

class _PARSER{
public: void init(){}
public: _PARSER(){}
public: std::vector <std::string>  splitstring(std::string work, char delim, int rep);
public: std::array <std::vector <std::string>,MAX_FILE_ENTRY> PARSEFILE(std::string filename, char DELIMITER);
public: std::array <std::vector <std::string>,MAX_FILE_ENTRY> PARSEFILE(std::string filename, char DELIMITER,string parsestr);
public: std::array <std::vector <std::string>,MAX_FILE_ENTRY> PARSEFILE1(std::string filename, char DELIMITER,string parsestr);
public: std::array<int,MAX_FILE_ENTRY> PARSEINT(std::array <std::vector <std::string>,MAX_FILE_ENTRY> parsedata, std::string parsestring, std::size_t nret);
public: std::array<double,MAX_FILE_ENTRY> PARSEDOUBLE(std::array <std::vector <std::string>,MAX_FILE_ENTRY> parsedata, std::string parsestring, std::size_t nret);
public: std::array<std::string,MAX_FILE_ENTRY> PARSESTRING(std::array <std::vector <std::string>,MAX_FILE_ENTRY> parsedata, std::string parsestring, std::size_t nret);

public :int PARSEINT(std::string filename, char delim, std::string searchstring);
public :double PARSEDOUBLE(std::string filename, char delim, std::string searchstring);
public :std::string PARSESTRING(std::string filename, char delim, std::string searchstring);
public :char PARSECHAR(std::string filename, char delim, std::string searchstring);
public :std::string remove_whitespace(std::string sstring);

};
#endif
