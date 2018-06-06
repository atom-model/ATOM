#ifndef CATMOSPHEREMODEL_H
#define CATMOSPHEREMODEL_H

#include <set>
#include <string>
#include <vector>

#include "tinyxml2.h"
#include "Array.h"
#include "Array_2D.h"

using namespace std;
using namespace tinyxml2;

class cAtmosphereModel {
public:
    const char *filename;

    cAtmosphereModel();
    ~cAtmosphereModel();

    // FUNCTIONS
    void LoadConfig(const char *filename);
    void Run();
    void RunTimeSlice(int time_slice);

    std::set<float>::const_iterator get_current_time() const{
        if(m_time_list.empty()){
            throw("The time list is empty. It is likely the model has not started yet.");
        }else{
            return m_current_time;
        }
    }

    std::set<float>::const_iterator get_previous_time() const{
        if(m_time_list.empty()){
            throw("The time list is empty. It is likely the model has not started yet.");
        }
        if(m_current_time != m_time_list.begin()){
            std::set<float>::const_iterator ret = m_current_time;
            ret--;
            return ret;
        }
        else{
            throw("The current time is the only time slice for now. There is no previous time yet.");
        }
    }

    bool is_first_time_slice() const{
        if(m_time_list.empty()){
            throw("The time list is empty. It is likely the model has not started yet.");
        }
        return (m_current_time == m_time_list.begin());
    }

    #include "AtmosphereParams.h.inc"

private:
    void SetDefaultConfig();

    std::set<float> m_time_list;
    std::set<float>::const_iterator m_current_time;
};

#endif
