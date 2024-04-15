#ifndef _SPMST_IO_H
#define _SPMST_IO_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex.h>

#if __cplusplus < 201703L
void inline read_string_param(std::istringstream &info){}
#endif

/**
 * @brief read parameters from stringstream
 * @example read_string_param(info,dt,dx,mystring)
 * 
 * @tparam Args you could read in many parameters
 * @param infile std::ifstream
 * @param args argments, arithmetic or string type 
 *
 */
template<typename T,typename ...Args>
void read_string_param(std::istringstream &info, T & first,Args& ...args)
{
    static_assert(std::is_arithmetic<T>::value ||
                    std::is_same<std::string,T>::value,
                    "please check input data type!\n");
    info >> first;
    #if __cplusplus < 201703L 
        read_string_param(info,args...);
    #else 
        if constexpr(sizeof...(args) > 0){
            read_string_param(info,args...);
        }
    #endif
}

/**
 * @brief read parameters from file
 * @example read_fmt_param(infile,dt,dx,mystring)
 * 
 * @tparam Args you could read in many parameters
 * @param infile std::ifstream
 * @param args argments, arithmetic or string type 
 *
 */
template<typename... Args>
void read_file_param(std::ifstream &infile,Args& ...args)
{
    std::string line;
    getline(infile,line);
    while(line.length() ==0 || line[0] == '#'){
        getline(infile,line);
    }
    std::istringstream info(line);
    read_string_param(info,args...);
    info.clear();
}

/**
 * @brief template function to read parameters by using regex
 * @example read_par_regex("NX",nx,"Par_file")
 * 
 * @tparam T datatype
 * @param filename filename
 * @param varname  variable name
 * @param var      reference to a variable
 */
template<typename T> int
read_par_regex(const std::string &varname,T &var,std::ifstream &infile)
{   
    // go to beginning of the file
    infile.clear();
    infile.seekg(0);

    // temporary vars
    std::string tempname = "^" + varname,dummy;
    std::istringstream info;
    int ierr = 1;

    // read line by line and find the match str
    std::string line;
    while(getline(infile,line)){
        if(line.length() == 0 || line[0] == '#') continue;;

        // match 
        regex_t reg;
        regmatch_t pmatch[1];
        regcomp(&reg,tempname.data(),REG_NEWLINE);
        int status = regexec(&reg,line.c_str(),1,pmatch,0);
        if(status == 0){
            info.str(line);
            info >> dummy; info >> dummy;
            info >> var;
            info.clear();
            ierr = 0;
            break;
        }
        regfree(&reg);
    }

    if(ierr == 1){
      printf("cannot find %s\n",varname.c_str());
    }
    
    return ierr;
};



/**
 * @brief template function to read string parameters by using regex
 * @example read_par_regex("NX",nx,"Par_file")
 * 
 * @tparam T datatype
 * @param filename filename
 * @param varname  variable name
 * @param var      reference to a variable
 */
template<> inline int 
read_par_regex<std::string>(const std::string &varname,std::string &var,
                            std::ifstream &infile)
{   
    // go to beginning of the file
    infile.clear();
    infile.seekg(0);

    // temporary vars
    std::string tempname = "^" + varname;
    int ierr = 1;

    // read line by line and find the match str
    std::string line;
    while(getline(infile,line)){
        if(line.length() == 0 || line[0] == '#') continue;;

        // match 
        regex_t reg;
        regmatch_t pmatch[1];
        regcomp(&reg,tempname.data(),REG_NEWLINE);
        int status = regexec(&reg,line.c_str(),1,pmatch,0);
        if(status == 0){
            size_t i = 0;
            for(; i < line.size(); i ++ ){
                if(line[i] == '=') break;
            }
            var = line.substr(i+1);
            ierr = 0;
            break;
        }
        regfree(&reg);
    }

    if(ierr == 1){
      printf("cannot find %s\n",varname.c_str());
    }
    
    return ierr;
}

void create_directory(const char *direcname);

#endif