#include "_constant.h"
#include "_declarations.h"
#include "_parser.h"
#include <algorithm>
#include <cctype>

std::vector <std::string> _PARSER :: splitstring(std::string work, char delim, int rep) {
   std::vector <std::string> flds;
   if (!flds.empty()) flds.clear();                                                         // empty vector if necessary
   std::string buf = "";
   std::size_t i = 0;
  
    while (i < work.length()) {

     if (work[i] != delim){
	   
       buf += work[i];
     }
     else if (rep == 1) {
       flds.push_back(this->remove_whitespace(buf));
       buf = "";
     } 
     else if (buf.length() > 0) {
       flds.push_back(this->remove_whitespace(buf));
       buf = "";
     }
     i++;
    }
    
    if (!buf.empty())
        flds.push_back(buf);
    return flds;
}

std::array <std::vector <std::string>,MAX_FILE_ENTRY> _PARSER :: PARSEFILE(std::string filename, char delim){
  std::string line;
  std::vector <std::string> flds;
  std:: ifstream fin;
  std::array <std::vector <std::string>, MAX_FILE_ENTRY> finaldata;
  std::size_t n=0;
  fin.open(filename,std::ifstream::in);
  if (fin.is_open()){
   while(getline(fin,line)) {
	   
    if (!line.empty()){
     flds = splitstring(line,delim,2);
      if (flds.size() > 1){
       for (std::size_t i=0; i<flds.size(); i++){
        finaldata[n].push_back(flds[i]);
       }
       n++;
      }
    }
   }
  }
  fin.close();
  return finaldata;
}

std::array <std::vector <std::string>,MAX_FILE_ENTRY> _PARSER :: PARSEFILE(std::string filename, char delim, std::string parsestr){
  std::string line;
  std::vector <std::string> flds;
  std:: ifstream fin;
  std::array <std::vector <std::string>, MAX_FILE_ENTRY> finaldata;
  
  //@ strings to parse for
  stringstream ss1,ss2;
  ss1<<"<"<<parsestr<<">";
  ss2<<"</"<<parsestr<<">";
  std::string startstr = ss1.str();
  std::string endstr = ss2.str();
  //@
  
  fin.open(filename,std::ifstream::in);
  bool proceedflag = false;
  std::size_t n=0;
  if (fin.is_open()){
    while(getline(fin,line)) {
	 cout<<proceedflag<<" "<<line<<endl;
     if (!line.empty()){
	 flds = this->splitstring(line,delim,2);
      if (flds.size() > 1 && proceedflag){
       for (std::size_t i=0; i<flds.size(); i++){
         finaldata[n].push_back(flds[i]);
       }
       cout<<endl;
       n++;
      }
      else if (flds.size() == 1) {
       if (flds[0].compare(startstr) == 0) proceedflag = true;
       if (flds[0].compare(endstr) == 0) proceedflag = false;
      }
     }
    }
  }
  fin.close();
  return finaldata;
}

std::array<int, MAX_FILE_ENTRY> _PARSER :: PARSEINT(std::array <std::vector <std::string>,MAX_FILE_ENTRY> parsedata, std::string parsestring, std::size_t nret){
	std::size_t ii = 0;
	std::array<int, MAX_FILE_ENTRY> retvalint;
	for (std::size_t i=0; i<parsedata.size(); i++){
	  if ( (!parsedata[i].empty()) && (parsedata[i][0].compare(parsestring) ==0)){
	    for (std::size_t j=1; j<nret+1; j++){
	     retvalint[ii] = std::atoi(parsedata[i][j].c_str());
	     ii += 1;
	    }
	  }
	}
	return retvalint;
}

std::array<double,MAX_FILE_ENTRY> _PARSER :: PARSEDOUBLE(std::array <std::vector <std::string>,MAX_FILE_ENTRY> parsedata, std::string parsestring, std::size_t nret){
	std::size_t ii = 0;
	std::array<double, MAX_FILE_ENTRY> retvaldoub;
	for (std::size_t i=0; i<parsedata.size(); i++){
	  if ((!parsedata[i].empty()) && (parsedata[i][0].compare(parsestring) ==0)){
	    for (std::size_t j=1; j<nret+1; j++){
	    	retvaldoub[ii] = std::atof(parsedata[i][j].c_str());
	    	ii += 1;
	    }
	  }
	}
	return retvaldoub;
}

int _PARSER :: PARSEINT(std::string filename, char delim, std::string searchstring){
  std::string line;
  std::vector <std::string> flds;
  std:: ifstream fin;
  int retval = 0;
  bool idenflag = false;
  fin.open(filename,std::ios::in);
  if (fin.is_open()){
   while(getline(fin,line) && !idenflag) {
    if (!line.empty()){
     flds = splitstring(line,delim,2);
      if ((flds.size() > 1) && (flds[0].compare(searchstring)==0)){              // if searchstring matches the read line
         retval = std::atoi(flds[1].c_str());
         idenflag = true;
        }
       }
      }
    }
  if (!idenflag) {
    cout<<"Could not find an entry for "<<searchstring<<endl;
    cout<<"Check parameter file "<<filename<<endl;
    cout<<"Exiting"<<endl;
    exit(1);
  }
  return retval;
}

double _PARSER :: PARSEDOUBLE(std::string filename, char delim, std::string searchstring){
  std::string line;
  std::vector <std::string> flds;
  std:: ifstream fin;
  double retval = 0;
  bool idenflag = false;
  fin.open(filename,std::ios::in);
  if (fin.is_open()){
   while(getline(fin,line) && !idenflag) {
    if (!line.empty()){
     flds = splitstring(line,delim,2);
      if ((flds.size() > 1) && (flds[0].compare(searchstring)==0)){              // if searchstring matches the read line
         retval = std::atof(flds[1].c_str());
         idenflag = true;
        }
       }
      }
    }
  if (!idenflag) {
    cout<<"Could not find an entry for "<<searchstring<<endl;
    cout<<"Check parameter file "<<filename<<endl;
    cout<<"Exiting"<<endl;
    exit(1);
  }
  return retval;
}

std::string _PARSER :: PARSESTRING(std::string filename, char delim, std::string searchstring){
  std::string line;
  std::vector <std::string> flds;
  std:: ifstream fin;
  std::string retval;
  bool idenflag = false;
  fin.open(filename,std::ios::in);
  if (fin.is_open()){
   while(getline(fin,line) && !idenflag) {
    if (!line.empty()){
     flds = splitstring(line,delim,2);
      if ((flds.size() > 1) && (flds[0].compare(searchstring)==0)){              // if searchstring matches the read line
         retval = this->remove_whitespace(flds[1].c_str());
         idenflag = true;
        }
       }
      }
    }
  if (!idenflag) {
    cout<<"Could not find an entry for "<<searchstring<<endl;
    cout<<"Check parameter file "<<filename<<endl;
    cout<<"Exiting"<<endl;
    exit(1);
  }
  
  return retval;
}

char _PARSER :: PARSECHAR(std::string filename, char delim, std::string searchstring){
  std::string line;
  std::vector <std::string> flds;
  std:: ifstream fin;
  char retval = '0';
  bool idenflag = false;
  fin.open(filename,std::ios::in);
  if (fin.is_open()){
   while(getline(fin,line) && !idenflag) {
    if (!line.empty()){
     flds = splitstring(line,delim,2);
      if ((flds.size() > 1) && (flds[0].compare(searchstring)==0)){              // if searchstring matches the read line
         const char* temp = flds[1].c_str();
         retval = temp[0];
         idenflag = true;
        }
       }
      }
    }
  if (!idenflag) {
    cout<<"Could not find an entry for "<<searchstring<<endl;
    cout<<"Check parameter file "<<filename<<endl;
    cout<<"Exiting"<<endl;
    exit(1);
  }
  return retval;
}

std::array<std::string,MAX_FILE_ENTRY> _PARSER :: PARSESTRING(std::array <std::vector <std::string>,MAX_FILE_ENTRY> parsedata, std::string parsestring, std::size_t nret){
	std::size_t ii = 0;
	std::array<std::string, MAX_FILE_ENTRY> retvalstr;
	for (std::size_t i=0; i<parsedata.size(); i++){
	  if ((!parsedata[i].empty()) && (parsedata[i][0].compare(parsestring) ==0)){
	    for (std::size_t j=1; j<nret+1; j++){
		    retvalstr[ii] = this->remove_whitespace(parsedata[i][j].c_str());
	    	    ii += 1;
	    }
	  }
	}
	return retvalstr;
}

std::string _PARSER :: remove_whitespace(std::string sstring){
  sstring.erase(std::remove(sstring.begin(), sstring.end(), ' '), sstring.end());
  return sstring;
}