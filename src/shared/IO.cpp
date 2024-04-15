#include <iostream>
#include <dirent.h>
#include <sys/stat.h>

/**
 * @brief Create a directory if it doesn not exist
 * 
 * @param dirname 
 */
void 
create_directory(const char *dirname)
{
    if(opendir(dirname) == NULL){
        mkdir(dirname,0777);
    } 
}